# Python module for reading a Torus parameters file
# D. Acreman, July 2018

class parameterList:
    """Class for list of TORUS parameters"""

    # The default file to read is parameters.dat
    def __init__(self, file='parameters.dat'):
        """Read the contents of a Torus parameters file."""
        
        # Create empty dictionaries and zero counters
        self.parameters = {}
        self.repeated_parameters = {}
        self.num_total_lines=0
        self.num_comment_lines=0
        self.num_blank_lines=0
        self.num_repeated_params=0
        self.filename=file
        
        print("Reading "+ self.filename)
        #  The with statement ensures the file is closed after being read
        with open(self.filename, "r") as file:
            for line in file:
                self.num_total_lines += 1
                # Strip out leading and trailing whitespace
                line_no_ws=line.strip()
                # Don't try to read the first element of an empty line
                if len(line_no_ws) > 0:
                # Discard comment lines
                    if line_no_ws[0] == '!':
                        self.num_comment_lines += 1
                    else:
                        # Split out part before any trailing comment
                        val=line_no_ws.split('!',maxsplit=1)[0]           
                        # Split out the first word of the string to use as the dictionary key
                        key = val.split()[0]
                        # If the key has already been read don't overwrite but record the information
                        if key in self.parameters:
                            print("WARNING: parameter "+key+" is repeated on line "+str(self.num_total_lines))
                            self.num_repeated_params += 1
                            if key in self.repeated_parameters:
                                self.repeated_parameters[key] += 1
                            else:
                                self.repeated_parameters[key] = 1
                        else:
                        # Add to the parameters dictionary using the first string as the key
                            self.parameters[key] = val.split()[1:]
                else:
                    self.num_blank_lines += 1

    def report_info(self):
        """ Report some information about the parameters in this file """
    
        print()
        print("========================================================")
        print("Filename:                         " + self.filename)
        print("Total number of lines read:       " + str(self.num_total_lines))
        print("Number of parameters:             " + str(len(self.parameters)))
        print("Instances of repeated parameters: " + str(self.num_repeated_params))
        print("Number of comment lines:          " + str(self.num_comment_lines))
        print("Numer of blank lines:             " + str(self.num_blank_lines))
        print("========================================================")
        print()

    def report_repeated_parameters(self):
        """Say which parameters were repeated and how many times"""

        print("Repeated parameters: ")
        for p in self.repeated_parameters:
            print(p, ' x ', str(self.repeated_parameters[p]))

    def report_integer_parameters(self):
        """Print a list of integer parameters and their values"""

        print("Integer parameters:")
        for p in self.parameters:
            try:
                print(p, int(self.parameters[p][0]))
            except ValueError:
                pass

# For when this module is run as a script
if __name__ == "__main__":
    p=parameterList()
    p.report_info()
    p.report_repeated_parameters()
    p.report_integer_parameters()
    
# End of file #####################################################################################################
