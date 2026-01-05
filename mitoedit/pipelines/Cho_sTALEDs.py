import logging
logger = logging.getLogger(__name__)
from .base_pipeline import BasePipeline

WINDOW_SIZE_MIN=14
WINDOW_SIZE_MAX=19
ADJACENT_BASES_LEFT_LEN=(16 + 15)
ADJACENT_BASES_RIGHT_LEN=(15 + 15)
MIN_LEFT_GAP=4
MIN_RIGHT_GAP=5
class ChosTALEDsPipeline(BasePipeline):
    """
    Cho sTALED Base Editing Pipeline (T→C and A→G edits)
    
    Implements the Cho et al. sTALED (split-TALED) base editing approach for mitochondrial DNA editing.
    This pipeline uses split transcription activator-like effector domains (sTALEDs) to target specific
    genomic locations for precise base editing, specifically converting T→C and A→G.
    
    Key Features:
    - Uses sTALED technology for targeted base editing
    - Supports T→C and A→G base conversions
    - Generates editing windows from both left and right sTALED positions
    - Handles multiple sequence contexts (CT, GT, TG, TC for T; AC, AG, CA, GA for A)
    - Optimized for mitochondrial DNA editing applications
    
    Unique Methods:
    - _generate_windows_from_right_sTALED(): Generates editing windows from right sTALED
    - _generate_windows_from_left_sTALED(): Generates editing windows from left sTALED
    """
    
    def __init__(self):
        super().__init__()
        self.pipeline_name = "Cho_G1397_sTALEDs"

    def _generate_windows_from_right_sTALED(self, circular_seq, pos, window_size):
        """Generate windows from right sTALED"""
        logger.debug("Generating windows from right sTALED.")
        windows = []
        for i in range(5, window_size - 4):
            start_index = pos - i
            end_index = pos + (window_size - i)
            window = circular_seq[start_index:end_index]
            windows.append(window)
        return windows

    def _generate_windows_from_left_sTALED(self, circular_seq, pos, window_size):
        """Generate windows from left sTALED"""
        logger.debug("Generating windows from left sTALED.")
        windows = []
        for i in range(5, window_size - 4):
            start_index = pos - (window_size - i)
            end_index = pos + i
            window = circular_seq[start_index:end_index]
            windows.append(window)
        return windows

    def process_mtDNA(self, mtDNA_seq, pos):
        """Main function which processes the DNA"""
        logger.info(f"Processing mtDNA sequence for position {pos}.")

        nospace_mtDNA = self._capitalize(self._remove_whitespace(mtDNA_seq))

        # list to store positions where 'T' is present in the mtDNA sequence in CT/GT/TG/TC contexts
        T_positions = self._find_dinucs(nospace_mtDNA, ['CT', 'GT', 'TG', 'TC'],'T')
        # list to store positions where 'A' is present in the mtDNA sequence in CA/GA/AG/AC contexts
        A_positions = self._find_dinucs(nospace_mtDNA, ['CA', 'GA', 'AG', 'AC'],'A')
        dummy = 0

        if pos in T_positions:
            ref, mut, all_windows,dum = 'T', 'C',[],[]
            circular_seq = nospace_mtDNA + nospace_mtDNA
            start_index = pos - ADJACENT_BASES_LEFT_LEN
            end_index = pos + ADJACENT_BASES_RIGHT_LEN
            adjacent_bases = circular_seq[start_index:end_index]
            T_positions_adj = self._find_dinucs(adjacent_bases, ['CT', 'GT', 'TG', 'TC'],'T')
            A_positions_adj = self._find_dinucs(adjacent_bases, ['CA', 'GA', 'AG', 'AC'],'A')
            #logger.info("The 60 adjacent bases to my target base are:", marked_adjacent)
            left_adjacent_bases = adjacent_bases[:30] #(base 0 to 29)
            right_adjacent_bases = adjacent_bases[31:] #(base 31 to 59)
            logger.info(f"The left and right adjacent bases are: {left_adjacent_bases} and {right_adjacent_bases}")
            FLAG=None
            #checking the position of the adjacent base --> if T present on either side --> then off-target
            if right_adjacent_bases[0] == 'T':
                dummy = 1
                dum.append(pos+1)
                FLAG=True
            if left_adjacent_bases[-1] =='T':
                dummy += 1
                dum.append(pos-1)
                FLAG=True
            if right_adjacent_bases[0] == 'A':
                FLAG=True
            if left_adjacent_bases[-1] =='A':
                FLAG=True
            for window_source in ["sTALED with AD on the left_TALE", "sTALED with AD on the right_TALE"]:
                for window_size in range(WINDOW_SIZE_MIN, WINDOW_SIZE_MAX): #for window sizes of 14-18bp long --> MAJOR ASSUMPTION!!
                    if window_source == "sTALED with AD on the left_TALE":
                        sTALED_left_T_windows = self._generate_windows_from_right_sTALED(circular_seq, pos, window_size) #generating the TC windows(14bp-18bp)
                        TALES=False
                        for rel_pos, window in enumerate(sTALED_left_T_windows, start=5):
                            # rel_pos is the position of target realive to the window
                            window_desc = f"Position {rel_pos} from the 5' end"
                            ws = f"{window_size}bp"
                            # marked_window = self._mark_bases(window, rel_pos, [(x + 3) for x in self._find_consecutive_GA_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_AC_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_AG_sequences(window[4:13])] + [(x + 3) for x in self._find_consecutive_CA_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_CT_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_GT_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_TG_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_TC_sequences(window[4:13])])
                            bystander_positions = [bystander_pos - ADJACENT_BASES_LEFT_LEN + rel_pos  for bystander_pos in (T_positions_adj + A_positions_adj) if bystander_pos != pos]
                            bystander_positions = [bystander_pos for bystander_pos in bystander_positions if bystander_pos > MIN_LEFT_GAP and bystander_pos <= window_size - MIN_RIGHT_GAP]
                            marked_window = self._mark_bases(window, rel_pos, bystander_positions)
                            print("left window",marked_window)
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos -rel_pos + 1
                            ac_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AC')
                            tc_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TC')
                            ag_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AG')
                            tg_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TG')
                            ca_positions = self._find_N_positions(window[3:12], start_position + 3, 'CA')
                            ct_positions = self._find_N_positions(window[3:12], start_position + 3, 'CT')
                            ga_positions = self._find_N_positions(window[3:12], start_position + 3, 'GA')
                            gt_positions = self._find_N_positions(window[3:12], start_position + 3, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return []
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))
    
                    elif window_source == "sTALED with AD on the right_TALE":
                        sTALED_right_T_windows = self._generate_windows_from_left_sTALED(circular_seq, pos, window_size) #generating the TC windows(14bp-18bp)
                        TALES=False
                        for rel_pos, window in enumerate(sTALED_right_T_windows, start=5):
                            window_desc = f"Position {rel_pos} from the 3' end"
                            ws = f"{window_size}bp"
                            #off_target_sites_set = count_TC_sequences(window[-13:-4]) | count_GA_sequences(window[-13:-4]) | count_AC_sequences(window[-13:-4]) | count_AG_sequences(window[-13:-4]) | count_CA_sequences(window[-13:-4]) | count_CT_sequences(window[-13:-4]) | count_GT_sequences(window[-13:-4]) | count_TG_sequences(window[-13:-4]) #finds the off-target sites in each window
                            #if dummy:
                            #    off_target_sites_set.add(dummy)
                            #off_target_sites = len(off_target_sites_set)
                            # marked_window = self._mark_bases(window, window_size - rel_pos + 1, [(x + window_size - 12) for x in self._find_consecutive_AG_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_AC_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GA_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CA_sequences(window[-13:-4])] + [(x + window_size - 12) for x in self._find_consecutive_TC_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_TG_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GT_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CT_sequences(window[-13:-4])])
                            bystander_positions = [bystander_pos - ADJACENT_BASES_LEFT_LEN + (window_size - rel_pos)  for bystander_pos in (T_positions_adj + A_positions_adj) if bystander_pos != pos]
                            bystander_positions = [bystander_pos for bystander_pos in bystander_positions if bystander_pos > MIN_LEFT_GAP and bystander_pos <= window_size - MIN_RIGHT_GAP]
                            marked_window = self._mark_bases(window, window_size - rel_pos, bystander_positions)
                            print("right window",marked_window)
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'T':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos - (window_size-rel_pos)
                            ac_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'AC')
                            tc_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'TC')
                            ag_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'AG')
                            tg_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'TG')
                            ca_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CA')
                            ct_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CT')
                            ga_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GA')
                            gt_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))

        elif pos in A_positions:
            logger.info("Base at position %d is in a editable context.", pos)
            ref, mut, all_windows, dum = 'A', 'G', [], []
            circular_seq = nospace_mtDNA + nospace_mtDNA
            start_index = pos - (16 + 15)
            end_index = pos + (15 + 15)
            adjacent_bases = circular_seq[start_index:end_index]
            T_positions_adj = self._find_dinucs(adjacent_bases, ['CT', 'GT', 'TG', 'TC'],'T')
            A_positions_adj = self._find_dinucs(adjacent_bases, ['CA', 'GA', 'AG', 'AC'],'A')
            left_adjacent_bases = adjacent_bases[:30] #(base 0 to 29)
            right_adjacent_bases = adjacent_bases[31:] #(base 31 to 59)
            logger.info(f"The left and right adjacent bases are: {left_adjacent_bases} and {right_adjacent_bases}")
            FLAG=None
            #checking the position of the adjacent base --> if T present on either side --> then off-target
            if right_adjacent_bases[0] == 'A':
                dummy = 1
                dum.append(pos+1)
                FLAG=True
            if left_adjacent_bases[-1] =='A':
                dummy += 1
                dum.append(pos-1)
                FLAG=True
            if right_adjacent_bases[0] == 'T':
                FLAG=True
            if left_adjacent_bases[-1] =='T':
                FLAG=True
            for window_source in ["sTALED with AD on the right_TALE", "sTALED with AD on the left_TALE"]: 
                for window_size in range(14, 19): #for window sizes of 14-18bp long
                    if window_source == "sTALED with AD on the right_TALE":
                        sTALED_left_A_windows = self._generate_windows_from_left_sTALED(circular_seq, pos, window_size) #generating the A windows(14bp-18bp)
                        TALES=False
                        for rel_pos, window in enumerate(sTALED_left_A_windows, start=5):
                            window_desc = f"Position {rel_pos} from the 3' end"
                            ws = f"{window_size}bp"
                            #off_target_sites = count_T_sequences(window[-13:-4]) + count_A_sequences(window[-13:-4])  #+dummy #finds the off-target sites in each window
                            # marked_window = self._mark_bases(window, window_size - rel_pos + 1, [(x + window_size - 12) for x in self._find_consecutive_AG_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_AC_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GA_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CA_sequences(window[-13:-4])] + [(x + window_size - 12) for x in self._find_consecutive_TC_sequences(window[-12:-3])] + [(x + window_size - 12) for x in self._find_consecutive_TG_sequences(window[-12:-3])] + [(x + window_size - 13) for x in self._find_consecutive_GT_sequences(window[-13:-4])] + [(x + window_size - 13) for x in self._find_consecutive_CT_sequences(window[-13:-4])])
                            bystander_positions = [bystander_pos - ADJACENT_BASES_LEFT_LEN + (window_size - rel_pos)  for bystander_pos in (T_positions_adj + A_positions_adj) if bystander_pos != pos]
                            bystander_positions = [bystander_pos for bystander_pos in bystander_positions if bystander_pos > MIN_LEFT_GAP and bystander_pos <= window_size - MIN_RIGHT_GAP]
                            marked_window = self._mark_bases(window, window_size - rel_pos, bystander_positions)
                            print("right window",marked_window)
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos - (window_size-rel_pos)
                            ac_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'AC')
                            tc_positions = self._find_N_positions(window[-12:-3], start_position +window_size- 12 - 1, 'TC')
                            ag_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'AG')
                            tg_positions = self._find_N_positions(window[-12:-3], start_position +window_size -12 - 1, 'TG')
                            ca_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CA')
                            ct_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'CT')
                            ga_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GA')
                            gt_positions = self._find_N_positions(window[-13:-4], start_position +window_size- 13, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))

                    elif window_source == "sTALED with AD on the left_TALE":
                        sTALED_right_A_windows = self._generate_windows_from_right_sTALED(circular_seq, pos, window_size) #generating the A windows(14bp-18bp)
                        TALES=False
                        for rel_pos, window in enumerate(sTALED_right_A_windows, start=5):
                            window_desc = f"Position {rel_pos} from the 5' end"
                            ws = f"{window_size}bp"

                            # marked_window = self._mark_bases(window, rel_pos, [(x + 3) for x in self._find_consecutive_GA_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_AC_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_AG_sequences(window[4:13])] + [(x + 3) for x in self._find_consecutive_CA_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_CT_sequences(window[3:12])] + [(x + 3) for x in self._find_consecutive_GT_sequences(window[3:12])] + [(x + 4) for x in self._find_consecutive_TG_sequences(window[4:13])] + [(x + 4) for x in self._find_consecutive_TC_sequences(window[4:13])])
                            bystander_positions = [bystander_pos - ADJACENT_BASES_LEFT_LEN + rel_pos  for bystander_pos in (T_positions_adj + A_positions_adj) if bystander_pos != pos]
                            bystander_positions = [bystander_pos for bystander_pos in bystander_positions if bystander_pos > MIN_LEFT_GAP and bystander_pos <= window_size - MIN_RIGHT_GAP]
                            marked_window = self._mark_bases(window, rel_pos, bystander_positions)
                            off_target_sites = marked_window.count('{') + dummy
                            int_pos_end = marked_window.find(']')
                            int_pos_ini = marked_window.find('[')
                            #to deal with T on either sides of the --> so it can be present in any context
                            if marked_window[int_pos_end + 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_end + 1)
                            elif marked_window[int_pos_ini - 1] == 'A':
                                final_window = self._mark_base_at_position(marked_window, int_pos_ini - 1)
                            else:
                                final_window = marked_window
                            start_position = pos -rel_pos + 1
                            ac_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AC')
                            tc_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TC')
                            ag_positions = self._find_N_positions(window[4:13], start_position + 3 , 'AG')
                            tg_positions = self._find_N_positions(window[4:13], start_position + 3 , 'TG')
                            ca_positions = self._find_N_positions(window[3:12], start_position + 3, 'CA')
                            ct_positions = self._find_N_positions(window[3:12], start_position + 3, 'CT')
                            ga_positions = self._find_N_positions(window[3:12], start_position + 3, 'GA')
                            gt_positions = self._find_N_positions(window[3:12], start_position + 3, 'GT')
                            ftc = tc_positions + tg_positions + ct_positions + gt_positions
                            ftg = ac_positions + ag_positions + ca_positions + ga_positions
                            ftg = [x for x in ftg if x != pos]
                            ftc = [x for x in ftc if x != pos]
                            combined_set = set(ftc + ftg)  # Combine into a set
                            sorted_combined = sorted(combined_set) if combined_set else []  # Sort if not empty, else return [0]
                            all_windows.append((self.pipeline_name, window_source, pos, ref, mut, ws, final_window, window_desc, off_target_sites, sorted_combined, TALES, FLAG))


        else:
            raise ValueError(f"Base at position {pos} is not in a editable context and cannot be edited by the {self.pipeline_name} pipeline.")

        return all_windows, adjacent_bases
