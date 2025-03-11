#!/usr/bin/env python3
import numpy as np
import subprocess as sp
from pathlib import Path

class Status:
    def __init__(self, system, username, submit_type, verbose=True):
        """ Initialize the Status class
        
        Parameters:
        -----------
        system : str
            The system to check
        username : str
            The username
        verbose : bool
            Increase output verbosity
        """
        self.system = system
        self.username = username
        self.verbose = verbose
        self.submit_type = submit_type
        if self.submit_type == "trial":
            self.submitfile = "submitted_jobs.log"
        elif self.submit_type == "equil":
            self.submitfile = "submitted_equil_jobs.log"
        elif self.submit_type == "prod":
            self.submitfile = "submitted_production_jobs.log"
        self.submitted = []
        return
    def check(self):
        """ Check the status of the jobs """
        self._obtain_submitted()
        self._obtain_running()
        self.check_status()
        #self._print_info(self)
    def _obtain_submitted(self):
        """ Obtain the list of submitted jobs """
            
        submitted = np.genfromtxt(f"{self.system}/{self.submitfile}", dtype=str)
        self.submitted = submitted
    def _obtain_running(self):
        """ Obtain the list of running jobs """
        running = sp.run(f'squeue -u {self.username} --format="%Z,%j,%t"', shell=True, stdout=sp.PIPE).stdout.decode('utf-8').split('\n')
        self.running = running
        self.print_list_sequential(running)
    def convert_jobname(self, jobname):
        """ Convert the jobname to the format used in the submitted jobs list 
        
        Parameters:
        -----------
        jobname : str
            The jobname to convert

        Returns:
        --------
        str
            The converted jobname
        """
        if self.submit_type == "trial":
            joblist = jobname.split('_')
            edge=joblist[0]
            sys=joblist[1]
            trial=joblist[4].split('.')[0]
            return f"{edge}_{sys}_t{trial}.slurm"
        elif self.submit_type == "equil":
            joblist = jobname.split('_')
            edge=joblist[1]
            sys=joblist[2]
            return f"{edge}_{sys}.slurm"
        elif self.submit_type == "prod":
            joblist = jobname.split('_')
            edge=joblist[1]
            sys=joblist[2]
            trial=joblist[5].split('.')[0]
            return f"{edge}_{sys}_t{trial}.slurm"
    def _check_pathname(self, pathname):
        """ Check if the pathname exists 
        
        Parameters:
        -----------
        pathname : str
            The pathname to check
        
        Returns:
        --------
        bool
            True if the pathname exists, False otherwise
            
        """
        if Path(pathname).exists():
            return True
        else:
            return False
    def _check_time_til_stepend(self, filename):
        """ Check the time remaining until the end of the step 
        
        Parameters:
        -----------
        filename : str
            The filename to check
        
        Returns:
        --------
        str
            The time remaining until the end of the step
        
        """
        with open(filename,'r') as f:
            for line in f:
                if "Estimated time remaining" in line:
                    return f"\n--> Step Time Left: {line.split(':')[-1].strip()}"
        return ""
    
    def check_status(self):
        """ Check the status of the jobs """
        system_path=Path(f"{self.system}/")
        count_complete = 0
        count_running = 0
        count_pd = 0
        for line in self.submitted:
            line = self.convert_jobname(line)
            found_match=False
            status = None
            for line2 in self.running:
                if line in line2 and system_path.name in line2:
                    found_match = True
                    status = line2.split(',')[2]

            if not found_match:
                if self.verbose: print(f"Job {line} has completed.")
                count_complete += 1
                if self.submit_type != "equil":
                    trialno = line.split('_')[2].split('.')[0]
                    if not self._check_pathname(system_path/f"unified/run/{line.split('_')[0]}/{line.split('_')[1]}/rem{trialno}.log"):
                        print(f"Job {line} has completed but the output file is missing.")
            if found_match:
                time_left = ""
                if status == "PD":
                    count_pd += 1
                    if self.verbose: print(f"Job {line} is pending.")
                elif status == "R" or status == "CD":
                    count_running += 1
                    try:
                        time_left = self._check_time_til_stepend(system_path/f"unified/run/{line.split('_')[0]}/{line.split('_')[1]}/mdinfo.000")
                    except:
                        time_left = ""
                    print(f"Job {line} is running, with status {status}. {time_left}")
        print("\n**********Summary:*********")
        print(f"----> Jobs completed: {count_complete}")
        print(f"----> Jobs running:   {count_running}")
        print(f"----> Jobs pending:   {count_pd}")
        print("***************************")
            
            
    def print_list_sequential(self, listname):
        """ Print the list sequentially
        
        Parameters:
        -----------
        listname : list
            The list to print
            
        """

        if self.verbose:
            for line in listname:
                print(line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", help="System to check. NOTE, system MUST be a directory in the current working directory.")
    parser.add_argument("--username", help="Username to check squeue info for.")
    parser.add_argument("--submit_type", default="trial", help="Submitted job type: Choose trial for equil_type=2, OR equil for equil_type=1, OR prod for equil_type=1")
    parser.add_argument("-v", "--verbose", help="Increase output verbosity, outputs all job details.", action="store_true")
    args = parser.parse_args()

    status = Status(args.system, args.username, args.submit_type, verbose=args.verbose)
    status.check()