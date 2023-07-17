import requests
from urllib.parse import urljoin
import time
import pprint
from typing import List, Tuple 
from abc import ABC, abstractmethod, abstractproperty
from typing import List, Dict

pp = pprint.PrettyPrinter(indent=4)

class Coster(ABC): 
    """ 
    Determines whether a molecule represented as a SMILES is 
    buyable. If it is, determines its cost 
    TODO: implement LookupCoster so user can input csv of inventory with cost=0
    or maybe allow ChemSpaceCoster to be updated with inventory list and 
    change those buyables to cost = 0
    """

    @abstractmethod
    def __call__(smis: List[str]) -> Tuple[Dict, List]: 
        """ For a list of smiles, returns a dictionary that maps the smiles to their cost
         and a list of the smiles that are buyable. """
    
    @abstractmethod
    def get_buyable_and_cost(smiles: str) -> Tuple[bool, float]:
        """ For a single smiles strings, outputs if the smiles is buyable and what its cost is"""

class NaiveCoster(Coster): 
    """ Determines if a molecule is buyable based on the number of certain elements. 
    Sets all costs to zero 
    """
    def __init__(self) -> None:
        return 
    
    def buyable(self, smiles: str) -> bool: 

        nc = smiles.count("c") + smiles.count("C")
        nn = smiles.count("N") + smiles.count("n")
        no = smiles.count("O") + smiles.count("o")
        
        if nc<=12 and nn<=3 and no<=5:
            return True        
        else: 
            return False
    
    def get_buyable_and_cost(self, smiles: str) -> Tuple[bool, float]: 

        buyable = self.buyable(smiles)

        if buyable: 
            return True, 0
        else: 
            return False, None 

    def __call__(self, smis: List[str]) -> Tuple[Dict, set]: 
        costs = {}
        buyables = set()

        for smiles in smis: 
            buyable = self.buyable(smiles)
            if buyable: 
                costs[smiles] = 0
                buyables.add(smiles)
            else: 
                costs[smiles] = None
                buyables.add(smiles)
        
        return costs, buyables 
                

class ChemSpaceCoster(Coster):
    """ Determines if a molecule is buyable and what its cost is (if buyable) 
    Adapted from PyChemSearch (Marco Stenta, Marta Pasquini @ Syngenta)"""
    def __init__(self, 
                 api_key, 
                 base_url='https://api.chem-space.com/', 
                 api_version='v3') -> None:

        self.api_key = api_key
        self.base_url = base_url
        self.api_version = api_version
        self.token_expiration_time = 3600
        self.get_token()
        self.categories_map = {
            # 'CSMS': 'Make-On-Demand Screening Compounds',
            # 'CSMB': 'Make-On-Demand Building Blocks',
            # 'CSCS': 'Custom Request', # only care about in stock
            'CSSB': 'In-Stock Building Blocks',
            'CSSS': 'In-stock Screening Compounds'
        }

        self.filters = {
            'purity': 90, 
            'shipsWithin': 10, 
            'min_packMg': 100, 
            'max_packMg': 10000,  
        } # TODO: make this user-set parameters

        return 

    def get_token(self):
        headers = {
            'Accept': 'application/json',
            'Authorization': 'Bearer {}'.format(self.api_key),
        }

        url = urljoin(urljoin(self.base_url, self.api_version), 'auth/token')
        r = requests.get(url, headers=headers, verify=False)
        response = r.json()

        if 'access_token' in response:
            token = response['access_token']
        else:
            print(response)
            token = None

        self.ref_time = time.time()  # set the reference time: the token can be used within
        self.token = token

    def single_search(self, smiles, max_results_count=1000,
                      categories='CSSB, CSSS'):

        # if the time difference between now and the reference time (set at token creation)
        # is larger than expiration time, a new token is requested
        delta_time = int(time.time() - self.ref_time)

        if delta_time > self.token_expiration_time:
            print("trigger token recalculation, time elapsed from token generation", delta_time)
            self.get_token()

        headers = {
            'Accept': 'application/json',
            'Authorization': 'Bearer {}'.format(self.token),
        }

        params = (
            ('count', max_results_count),
            ('categories', categories),
        )

        files = {
            'SMILES': (None, smiles),
        }

        url = urljoin(self.base_url, '{}/search/{}'.format(self.api_version, 'exact'))

        response = requests.post(url, headers=headers, params=params, files=files,  verify=False)
        processed_response = self.process_response(response)
        
        return processed_response

    def cost_from_response(self, response):
        offers = [offer for item in response['content']['items'] for offer in item['offers']]
        
        if len(offers) == 0: 
            return None 
        
        filtered_offers = []

        for offer in offers: 
            if self.offer_filter(offer): 
                [filtered_offers.append(price) for price in offer['prices']]

        prices = [offer['priceUsd']*1000/offer['packMg'] for offer in filtered_offers if self.price_filter(offer)] # USD/g

        return min(prices)

    def offer_filter(self, offer):
        if 'shipsWithin' in self.filters and int(offer['shipsWithin']) > self.filters['shipsWithin']: 
            return False
        if 'purity' in self.filters and offer['purity'] < self.filters['purity']: 
            return False 
        return True 
    
    def price_filter(self, price): 
        if 'min_packMg' in self.filters and price['packMg'] < self.filters['min_packMg']: 
            return False 
        if 'max_packMg' in self.filters and price['packMg'] > self.filters['max_packMg']: 
            return False
        if 'max_rawCost' in self.filters and price['priceUsd'] > self.filters['max_rawCost']:
            return False 
        return True 

    def process_response(self, response): 
        if response.status_code == 200:
            content = response.json()
        else:
            content = None
            print('msg from process_response():', 'status_code', response.status_code, 'reason', response.reason)

        return {'content': content, 'status_code': response.status_code, 'reason': response.reason}
    
    def get_buyable_and_cost(smiles: str) -> Tuple[bool, float]:
        
        response = self.single_search(smiles)
        cost = self.cost_from_response(response)
        if cost is None: 
            return False, None
        else: 
            return True, cost

    def __call__(self, smis: List[str]): 

        costs = {}
        buyables = set()

        for smiles in smis: 
            response = self.single_search(smiles)
            cost = self.cost_from_response(response)
            costs[smiles] = cost
            if cost is not None: 
                buyables.add(smiles)

        
        return costs, buyables 