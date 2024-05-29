""" This module contains the Selector ABC and its implementations. A selector orchestrates SPARROW's downselection by 
collecting all route graph information, formulating an optimization problem, solving it, and outputting results. """

from sparrow.selector.base import Selector