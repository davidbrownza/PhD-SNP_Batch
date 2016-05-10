#!/usr/bin/env python

import os, sys, argparse, subprocess, datetime, time


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False 


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False 


class PhDSNPResult:
    
    def __init__(self, name, prediction, reliability_index):
        self.name = name
        self.prediction = prediction
        self.reliability_index = reliability_index
    
    def __str__(self):
        return "%s\t%s  \t%s" % (self.name, self.prediction, self.reliability_index)
    
    
class PhDSNP:
    
    def __init__(self, mode, seq=None, pos=None, new_res=None):
        self.mode = mode
        self.seq = seq
        self.pos = pos
        self.new_res = new_res
        self.results = None


    def validate(self):
        assert self.mode in ["-seq", "-seqp", "-seqb"], "invalid mode"
        assert os.path.exists(self.seq), "sequence file does not exist"
        assert is_integer(self.pos), "position must be an integer"
        assert self.new_res.isalpha(), "new residue must be a letter character"
    
    
    def compile_command(self):
        phd_exe = os.path.join("${PHDHOME}", "PhD-SNP.py")
        
        command = "python -O %s %s %s %s %s" % (
                phd_exe, self.mode, self.seq, self.pos, self.new_res
            )
            
        return command
    
    
    def parse_result(self, out, err):
        sys.stderr.write(err)
        
        scores = []
        started = False
        
        lines = out.split("\n")
        for line in lines:
            if not started and (line.startswith("      Profile") or line.startswith("      Sequence")):
                started = True
            
            elif started and len(line) > 0:
                pos = line[7:14].strip()
                wild_type = line[18].strip()
                new_type = line[23].strip()
                    
                name = "%s%s%s" % (wild_type, pos, new_type)
                
                prediction = line[31:39].strip()
                reliability_index = line[41:43].strip()
                
                scores.append(PhDSNPResult(name, prediction, reliability_index))
            
            elif started:
                break
        
        self.results = scores
        return self.results


    def submit(self):
        cmd = self.compile_command()
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return process


    @staticmethod
    def load_from_file(path):
        params = []
        
        with open(path, 'r') as f:
            count = 0
            for line in f:
                try:
                    count += 1
                    
                    arr = line.split("|")
                    assert len(arr) >= 4, "not enough parameters were provided"
                    
                    phd = PhDSNP(mode=arr[0], seq=arr[1], pos=arr[2], new_res=arr[3].strip())
                    phd.validate()
                        
                    params.append(phd)
                        
                except Exception, ex:
                    sys.stderr.write("Line %d: %s\n" % (count, str(ex)))
                
        return params


def kill_process(proc):
    try:
        proc.kill()
    except OSError:
        # proc is already dead
        pass


def run(params, num_processes=1):
    processes = []
    scores = []
    
    submitted = 0
    completed = 0
    while completed < len(params):
        
        if len(processes) < num_processes:
            # spawn a new process
            phd = params[submitted]
            proc = phd.submit()
            
            processes.append((phd, proc))
            
            submitted += 1
        
        else:
            # monitor running processes
            while len(processes) >= num_processes or submitted == len(params):
                
                for process in processes:
                    proc = process[1]
                    
                    # if process completed
                    if proc.poll() is not None:
                        out, err = proc.communicate()
                        
                        phd = process[0]
                        scores += phd.parse_result(out, err)
                        
                        # kill process
                        kill_process(proc)
                        
                        completed += 1
                        print "Completed %d/%d" % (completed, len(params))
                
                processes = [process for process in processes if process[1].returncode is None]
                if len(processes) == 0:
                    break
    
    return scores
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Submit multiple mutations to I-Mutant 2.0")
    
    parser.add_argument('--processes', type=int, help='Number of processes to run in parallel', default=1)
    parser.add_argument('--input', help='A batch input file consisting of the parameters to submit to the I-Mutant 2.0 script')
    parser.add_argument('--output', help='Output file name')
    
    args = parser.parse_args(sys.argv[1:])
    
    params = PhDSNP.load_from_file(args.input)
    num_processes = min(len(params), args.processes) 
    print "Running PhD-SNP on %d cores..." % num_processes

    start_time  = datetime.datetime.now()
    scores = run(params, num_processes)
    end_time = datetime.datetime.now()
    
    time_taken = end_time - start_time
    print "Time taken:" + str(time_taken)
    
    with open(args.output, "w") as f:
        f.write("#Name\tPrediction \tRI")
        
        for result in scores:
            f.write("\n%s" % result)
    
