import re
import shutil
import os

class CheckNormalT():
    '''
    Check if every file in HF and CI folders is complete

    Methods
    -------
    check(move=True)
    Check files

    move (`bool`):
        If True, move files inside folder 'Incomplete'
    '''
    def __init__(self):
        self.base_path = os.getcwd()
        self.NORMAL = 'Normal termination'
        
        # Valuable re variables
        self.normal_re = re.compile(self.NORMAL)

        sets = [os.path.join(self.base_path, set) for set in os.listdir(self.base_path) if 'SET' in set]
        self.sets = sorted(sets)

    def check(self, move=True):
        '''
        Check files

        Parameters
        ----------
        move (`bool`):
            If True, move files inside folder 'Incomplete'
        '''
        print('Checking all files before run calculations...\n')
        
        for path in self.sets:
            HFt = []
            CIt = []
            
            for Dir, SubDir, Files in os.walk(path):
                if '80-e' in Dir or 'incomplete' in Dir:
                    continue
                elif 'HF' in Dir:
                    HF_p = Dir
                    for file in Files:
                        if 'log' in file:
                            HFt.append(os.path.join(Dir, file))
                elif 'CI' in Dir:
                    CI_p = Dir
                    for file in Files:
                        if 'log' in file:
                            CIt.append(os.path.join(Dir, file))
                
            for file in HFt:

                    # Read all lines
                    f = open(file, 'r')
                    lines = f.readlines()
                    f.close()

                    search = str(self.normal_re.search(str(lines[len(lines)-1])))

                    if search == 'None':
                        file_name = re.search(r'[/\\]([A-Z0-9]+[_a-z]*.log)', file).group(1)

                        destination_p = os.path.join(HF_p, 'incomplete_calculations')

                        try:
                            os.mkdir(destination_p)
                        except Exception:
                            pass

                        print ('Incomplete gaussian calculation for', file_name)

                        source = file
                        destin = os.path.join(destination_p, file_name)

                        if move:
                            shutil.move(source, destin)

            for file in CIt:

                    # Read all lines
                    f = open(file, 'r')
                    lines = f.readlines()
                    f.close()

                    search = str(self.normal_re.search(str(lines[len(lines)-1])))

                    if search == 'None':
                        file_name = re.search(r'[/\\]([A-Z0-9]+[_a-z]*.log)', file).group(1)

                        destination_p = os.path.join(CI_p, 'incomplete_calculations')

                        try:    
                            os.mkdir(destination_p)
                        except Exception:
                            pass

                        print('Incomplete gaussian calculation for', file_name)

                        source = file
                        destin = os.path.join(destination_p, file_name)
                        
                        if move:
                            shutil.move(source, destin)