import subprocess
import os
import isoform2


def run_hmm(hmm, fasta):
    ''' cmd to run hmm2'''
    dir     = os.path.dirname(hmm)
    models  = os.path.join(dir, "models")
    
    input = [
        os.path.join(models, "don.pwm"),
        os.path.join(models, "acc.pwm"),
        os.path.join(models, "exon.mm"),
        os.path.join(models, "intron.mm"),
        os.path.join(models, "exon.len"),
        os.path.join(models, "intron.len")
    ]
    
    cmd = [hmm, fasta] + input
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    
    return result.stdout

def run_geniso2(geniso, fasta, model, hints=False):
    ''' cmd to run geniso '''
    
    cmd   = [geniso, fasta, model]
    if hints:
        cmd.append('--hmm')
        cmd.append(str(hints))
    try:
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print(f"STDERR: {e.stderr}")
        # Continue execution without failing by returning the error message
        return f"Error in geniso2: {e.stderr}"

def parse(output):
    ''' Parse HMM output to extract donor and acceptor site probabilities '''
    
    dons = []
    accs = []
    switch = 0
    
    lines = output.strip().split('\n')
    for line in lines:
        # this switch logic is odd, lit can turn into anything better
        if      line == "DONS": 
            switch = 0
            continue
        elif    line == "ACCS":
            switch = 1
            continue
        
        line = line.strip().split('\t')
        pos = int  (line[0])
        val = float(line[1])

        if   switch == 0: dons.append( (pos, val) )
        elif switch == 1: accs.append( (pos, val) )
    
    return dons, accs

def gapstats(sites):
    ''' Gap statistic to find significant splice sites '''
    
    sites = sorted(sites, key=lambda x: x[1], reverse=True)
    max = 0
    idx = 0
    
    for i in range( len(sites) - 1 ):
        val  = sites[i][1]
        val2 = sites[i+1][1]
        gap  = val / val2

        if gap > max:
            max = gap 
            idx = i
    
    return sorted( [site[0] for site in sites[:idx+1]] )

def smoothed_gapstats(sites, k:int = 5):
    ''' smoothed gap statistics '''
    
    sites= sorted(sites, key=lambda x: x[1], reverse=True)
    max = 0
    idx = 0

    for i in range( len(sites) - k ):
        avg  = sum(site[1] for site in sites[i:i+k])     / k
        avg2 = sum(site[1] for site in sites[i+1:i+k+1]) / k
        gap  = avg/avg2 

        if gap > max:
            max = gap
            idx = i + k // 2
    
    return sorted( [site[0] for site in sites[:idx+1]] )

def percentile(sites, percentile:int = 25):
    ''' filter splice sites by percentile '''
    
    sites  = sorted(sites, key=lambda x:[1], reverse=True)
    cutoff = sites[int ( len(sites) * (1 - percentile/100) )][1]
    sites  = [site for site in sites if site[1] >= cutoff]
    return sorted([site[0] for site in sites])

def countiso(dons, accs, min_intron, min_exon, limit=False):
    ''' count all possible isoform '''
    
    count   = 0
    def all_possible(dpos, apos, ipos, old):
        ''' recursion '''
        nonlocal count
        
        if limit:
            if count > limit: return
            
        test  = ipos == 0
        start = dpos      if test else apos
        end   = len(dons) if test else len(accs)
        
        for i in range(start, end):
            new = dons[i] if test else accs[i]
            if old != 0 and test:      
                if new - old - 1 < min_exon: continue
            else:                       
                if new - old + 1 < min_intron: continue
            if not test: count += 1
            if test:    all_possible(i + 1, apos, 1, dons[i])
            else:       all_possible(dpos, i + 1, 0, accs[i])
                
    all_possible(0, 0, 0, 0)
    return count