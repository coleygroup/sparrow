from abc import ABC, abstractmethod
from typing import Dict
import subprocess
import pandas as pd
from pathlib import Path
import time 

class RxnClass(ABC):
    def __init__(self):
        self.status_log = None
        self.no_class_num = 0

    @abstractmethod
    def get_rxn_class(rxn: str) -> str:
        """Maps single reaction to its class, if it has one. If not, returns string 0."""

    @abstractmethod
    def get_rxn_classes(rxn_file: Path) -> Dict:
        """Converts a list of reactions into a dictionary of keys, corresponding to reaction classes, which map to lists of reactions"""

class NameRxnClass(RxnClass):
    def __init__(self, 
                 dir: str = None,
                 tmp: str = './tmp/'):
        self.dir = dir
        self.tmp = Path(tmp)
        self.tmp.mkdir(parents=True, exist_ok=True)
        super().__init__()

    def get_rxn_class(self, rxn, attempt=0):
        # need to test this 
        if attempt > 5:
            print(f"Classifying {rxn} failed")
            self.no_class_num += 1
            return f"Unclassified_{self.no_class_num}"
        
        if rxn[0].startswith('>'): 
            self.no_class_num += 1
            return f"Unclassified_{self.no_class_num}"
        
        tmp_in = self.tmp / "rxn.smi"

        with open(tmp_in, "w") as f: 
            f.write(rxn)

        cmd_str = " ".join([self.dir + "/namerxn", "-nomap", str(tmp_in)])
        output = subprocess.check_output(cmd_str, shell=True).decode("utf-8").split('\n')[0]
        class_num = self.output_to_classnum(output)
        # class_num = ""
        # with open(tmp_out, "r") as f:
        #     output = f.readlines()
        # # output = open("test.smi.out", 'r')

        # try: 
        #     class_num = output[0].strip().split()[1]
        # except:
        #     return self.get_rxn_class(rxn, attempt + 1)
            
        return class_num

    def get_rxn_classes(self, rxns):

        tmp_in = self.tmp / "rxns.smi"
        with open(tmp_in, "w") as f: 
            f.write('\n'.join(rxns))
        # tmp_out = self.tmp / "rxn.smi.out"
        # cmd_str = " ".join(["../" + self.dir + "/namerxn", "-completer", "-addrxnname", "-osmi", str(rxn_file), "-o", str(tmp_out)])
        # subprocess.Popen(cmd_str, shell=True)

        cmd_str = " ".join([self.dir + "/namerxn", "-nomap", "-addrxnname", "-osmi", str(tmp_in)])
        output = subprocess.check_output(cmd_str, shell=True).decode("utf-8").split('\n')

        classes = [self.output_to_classnum(line) for line in output[:len(rxns)]]
        
        return classes 
    
        # with open(tmp_out, "r") as f: 
        #     output = f.readlines()

        # for line in output[:rxns]: 
        #     class_num = ""
        #     line = line.strip().split()
        #     if len(line) >= 2:
        #         class_num = line[1]
        #     if class_num == "" or class_num == None:
        #         self.no_class_num += 1
        #         class_num = "Unclassified by NameRxn" + str(self.no_class_num)

        #     if class_num not in classes.keys():
        #         classes[class_num] = []
        #     if line != []:
        #         classes[class_num].append(line[0])
        # return classes
    
    def output_to_classnum(self, line): 
        if line == '': 
            self.no_class_num += 1
            class_num = f"Unclassified_{self.no_class_num}"
        
        try: 
            class_num = line.split(' ')[1]
        except: 
            self.no_class_num += 1
            class_num = f"Unclassified_{self.no_class_num}"

        if class_num.startswith('0'): 
            self.no_class_num += 1
            class_num = f"Unclassified_{self.no_class_num}"    

        return class_num 
    

class LookupClass(RxnClass):
    def __init__(self, 
                 csv_path: str = None):
        # TODO: test the csv -> dict conversion
        super().__init__()
        
        df = pd.read_csv(csv_path)
        self.dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
        self.no_class_num = self.update_no_class_num()

    def get_rxn_class(self, rxn):
        if self.dict: 
            return self.dict.get(rxn, "None")
        self.no_class_num += 1
        return f"Unclassified_{self.no_class_num}"    

    def get_rxn_classes(self, rxns):
        classes = [self.get_rxn_class(rxn) for rxn in rxns]
        return classes

    def update_no_class_num(self):  
        for value in self.dict.values():
            if value.startswith("Unclassified_"): 
                try: 
                    unclass_num = float(value.strip("Unclassified_"))
                except: 
                    continue 
                self.no_class_num = max(self.no_class_num, unclass_num)
        
        return int(self.no_class_num)