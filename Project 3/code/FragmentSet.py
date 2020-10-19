import os
import re


class FragmentSet(object):
    def __init__(self, fragfile, rmsdfile):
        """
        This class contains the fragment library for the input protein. It must do the following:
        - Read in fragment file and parse fragments at each position. Fragment files are of the form <protein>_<frag_len>mers.frag
        - Read in RMSD file containing pre-calculated RMSD to native structure for each fragment at each position.
        - Based on fragments and their corresponding RMSDs, rank fragments at each position by RMSD
        """
        try:
            with open(fragfile) as f:
                self.fragfile = self.__preprocessFragFile__(f.read())
        except FileNotFoundError:
            print("Fragment file not found.")
        try:
            with open(rmsdfile) as f:
                self.rmsdfile = self.__preprocessRmsdFile__(f.read())
        except FileNotFoundError:
            print("RMSD file not found.")
        return

    @staticmethod
    def __preprocessRmsdFile__(raw_rmsd_file: str) -> dict:
        lines = raw_rmsd_file.split('\n')
        rmsd_dict = {}
        pos = 1
        sub_dict = {}
        for line in lines:
            if line:
                elems = line.split('\t')
                new_pos = int(elems[0])
                # if new position
                if new_pos != pos:
                    # add dict to outer dict
                    rmsd_dict[pos] = sub_dict
                    # re init things
                    pos = new_pos
                    sub_dict = {}
                # add this line
                sub_dict[int(elems[1])] = float(elems[2])
        return rmsd_dict

    def __preprocessFragFile__(self, raw_frag_file: str) -> dict:
        """
        Process fragment file to be of the format:
        dict of list(position) of tuple of int(index of fragments) and list(tuple of phi and psi)
        e.g. {pos: (index, [(180, 180), ...])}
        :param raw_frag_file: file content as a string
        :return: formatted fragment file
        """
        lines = raw_frag_file.split('\n')
        frag_dict = {}
        row_num = 0
        index = 0
        cur_lst = []
        for row in lines:
            if row.startswith(" position:"):
                row_num = self.__findPosition__(row)
                frag_dict[row_num] = []
                cur_lst = []
            elif row:
                index = int(row[92:94])
                cur_lst.append((float(row[19:27]), float(row[28:36])))
            else:
                if cur_lst:
                    frag_dict[row_num].append((index, cur_lst))
                cur_lst = []
        return frag_dict

    @staticmethod
    def __findPosition__(row: str) -> int:
        position = re.search('position:            ([0-9]+) neighbors', row)
        try:
            result = int(position.group(1))
            return result
        except AttributeError:
            print("Cannot find position.")
            return 0

    def get_lowRMS_fragments(self, pos, N):
        """
        Returns the top-ranked fragments by RMSD at a defined position in the chain
        --------
        Params
            - pos (int): fragment position in chain (1-indexed)
            - N (int): number of fragments to return
        Returns
            - lowRMS_fragments (list): top N fragments at pos by RMSD. This should be a list of lists of (phi, psi) tuples.
              For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]
        """
        # get indexes of least rmsd fragments
        rmsd_dict = self.rmsdfile[pos]
        selected_index = sorted(rmsd_dict, key=rmsd_dict.get)[:N]
        # return those fragments
        return [t[1] for t in self.fragfile[pos] if t[0] in set(selected_index)]


def debug():
    from pathlib import Path
    cur_path = os.getcwd()
    project_dir = Path(cur_path).parent
    data_dir = os.path.join(cur_path, "starter_data")
    frag = "helix_9mers.frag"
    rmsd = "helix_9mers.rmsd"
    fs = FragmentSet(os.path.join(data_dir, frag), os.path.join(data_dir, rmsd))
    fs.get_lowRMS_fragments(3, 3)


# if __name__ == '__main__':
#     debug()
