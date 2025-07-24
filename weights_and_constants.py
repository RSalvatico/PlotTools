import numpy as np

class weights_and_constants:
    """
    This class contains the weights and constants used in the analysis.
    """

    def __init__(self):
        self.weights_0p6ttWcb_and_0p05ttLF = {
            "ttLF": 0.537,
            "ttcc": 0.09,
            "ttcj": 0.116,
            "ttbb": 0.071,
            "ttbj": 0.156
        }
        self.weights_0p6ttWcb_and_0p1ttLF = {
            "ttLF": 0.63,
            "ttcc": 0.09,#should be 0.01 technically
            "ttcj": 0.12,
            "ttbb": 0.04,
            "ttbj": 0.10
        }

        # Define event classification selection and binning

        #################
        evtClassification_weights = self.weights_0p6ttWcb_and_0p1ttLF # Change here the set of weights to use
        #################

        eventClassificationBaseSelection = "score_tt_Wcb > 0.6 && score_ttLF < 0.1"
        SR_selection = "score_tt_Wcb > 0.85"
        CR_selection = "score_tt_Wcb < 0.85"
        self.adhoc_selection = {
            "score_tt_Wcb" : f"{eventClassificationBaseSelection} && {SR_selection}",
            "fscore_ttbb"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttcj']} * score_ttcj && {evtClassification_weights['ttbb']} * score_ttbb > {evtClassification_weights['ttLF']} * score_ttLF",
            "fscore_ttbj"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttcj']} * score_ttcj && {evtClassification_weights['ttbj']} * score_ttbj > {evtClassification_weights['ttLF']} * score_ttLF",
            "fscore_ttcc"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttcj']} * score_ttcj && {evtClassification_weights['ttcc']} * score_ttcc > {evtClassification_weights['ttLF']} * score_ttLF",
            "fscore_ttcj"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttcj']} * score_ttcj > {evtClassification_weights['ttLF']} * score_ttLF",
            "fscore_ttLF"  : f"{eventClassificationBaseSelection} && {CR_selection} && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttbb']} * score_ttbb && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttbj']} * score_ttbj && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttcc']} * score_ttcc && {evtClassification_weights['ttLF']} * score_ttLF > {evtClassification_weights['ttcj']} * score_ttcj"
        }
        self.adhoc_binning = {
            "score_tt_Wcb" : np.array([0.,0.9,1.]),
            "fscore_ttbb"  : np.array([0.,0.7,1.]),
            "fscore_ttbj"  : np.array([0.,0.45,1.]),
            "fscore_ttcc"  : np.array([0.,0.45,1.]),
            "fscore_ttcj"  : np.array([0.,0.35,1.]),
            "fscore_ttLF"  : np.array([0.,0.1,1.]),
        }

_wc_instance = weights_and_constants()
adhoc_selection = _wc_instance.adhoc_selection
adhoc_binning = _wc_instance.adhoc_binning