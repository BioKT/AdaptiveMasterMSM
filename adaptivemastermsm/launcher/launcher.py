"""
This file is part of the AdaptiveMasterMSM package.

"""
#!/usr/bin/env python

import sys, os
import subprocess
import shlex
import multiprocessing as mp

# AdaptiveMasterMSM
from adaptivemastermsm.system import system

class Launcher(object):
    """
    Call GROMACS and process output to pass it to class analyzer
    """

    def __init__(self, md_step, forcefield, water, \
            filepdb, wet_xtc_file,  dry_xtc_file=None, last_wet_snapshot=None):
        """
        Read input from system, run MD, and process output
        for the class analyzer

        Args:
            sysfile (str): File containing system information?
            md_step (str): 'Equilibration' or 'Production'
            water: water model for GROMACS (tip3p, ...)
            filepdb (str): PDB file defining the system
            wet_xtc_file (str)
            dry_xtc_file (str)
            last_wet_snapshot (str)

        Return?:
            trajfile (str): File containing processed trajectories
        """
        self.pdb = filepdb
        self.dry_xtc_file = dry_xtc_file
        self.last_wet_snapshot = last_wet_snapshot
        # create an instance of system class to define the system and parameters
        params = system.System(water, md_step, 1)

        if params.md_step == 'Production':
            self.wet_xtc_file = wet_xtc_file           
            n_short_runs = 2 # this will be managed by controller
            mdp_template = params.filemdp

            # ndx = ??? use this to set positions of ligands in future work
            for i in range(n_short_runs):
                mdp = self.gen_mdp(mdp_template, i)
                # different resampling positions:
                gro = "data/gro/%s.gro" % i
                if not os.path.exists(gro[:gro.rfind("/")]):
                    os.makedirs(gro[:gro.rfind("/")])
                cmd = "cp %s data/gro/%s.gro" % (self.pdb, i)
                os.system(cmd)
                top = "data/top/%s.top" % i
                if not os.path.exists(top[:top.rfind("/")]):
                    os.makedirs(top[:top.rfind("/")])
                itp = "data/itp/%s.itp" % i
                if not os.path.exists(itp[:itp.rfind("/")]):
                    os.makedirs(itp[:itp.rfind("/")])
                cmd = 'gmx pdb2gmx -f %s -o processed_%s -ff %s -water %s' % \
                    (gro, self.pdb, forcefield, params.water)
                print(" running: ",cmd)
                p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                output, error = p.communicate()
                #print(output, error)
                cmd = "mv processed_%s %s & mv topol.top %s & mv posre.itp %s" % \
                        (self.pdb, gro, top, itp)
                os.system(cmd)
                out = "data/tpr/%s.tpr" % i
                output, error = self.gen_tpr(gro, top, mdp, out, i)

            # set multiprocessing options
            n_threads = mp.cpu_count()
            pool = mp.Pool(processes=n_threads)
            # run simulations
            gmxinput = [[out, gro] for i in range(n_short_runs)]
            results = []
            for x in gmxinput:
                results.append(pool.apply_async(self.gromacs_worker, [x]))
            # close the pool and wait for each running task to complete
            pool.close()
            pool.join()
            for result in results:
                out, err = result.get()
                print("out: {} err: {}".format(out, err))

        elif params.md_step == 'Equilibration':
            # first equilibrate
            cmd = 'gmx pdb2gmx -f %s -o processed_%s -ff %s -water %s -ignh' % \
                    (self.pdb, self.pdb, forcefield, params.water)
            print(" running: ",cmd)
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            #print(out, err)
            out, err = params.build_box(self.pdb, 1)
            #print(out, err)
            user_wet_xtc_file = wet_xtc_file
            self.wet_xtc_file = 'equilibration.xtc'
            self.run_md(4, params)

            # then production
            # p.md_step = 'Production'  #ionix, ojo aqui, por eso i+1=2
            p_prod = system.System(water, 'Production', 2)
            self.wet_xtc_file = user_wet_xtc_file
            self.run_md(4, p_prod)

        self.clean_working_directory()

        return

    def gen_tpr(self, gro, top, mdp, tpr, i):
        """
        Function for generating tpr files

        """
        if not os.path.exists(tpr[:tpr.rfind("/")]):
            os.makedirs(tpr[:tpr.rfind("/")])
        cmd = "gmx grompp -c %s -p %s -f %s -o %s -maxwarn 1"\
                %(gro, top, mdp, tpr)
        print(" running: ",cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        # print(out, err)
        cmd = "mv mdout.mdp data/mdp/%s_out.mdp" % i
        os.system(cmd)
        return out, err
    
    def gen_mdp(self, filemdp, i):
        """ 
        Generate Gromacs parameter input file (mdp)
        
        """
        filemdpout = "data/mdp/%s.mdp" % i
        try:
            out = open(filemdpout, "w")
        except IOError:
            os.makedirs(filemdpout[:filemdpout.rfind("/")])
            out = open(filemdpout, "w")
        rawmdp = open(filemdp).readlines()
        for i in rawmdp:
            out.write("%s"%i)
        out.close()
        return filemdpout
    
    def gromacs_worker(self, x):
        """
        Worker function for running Gromacs jobs
        
        md_step = x[0]
        tpr = x[1]
        mdp = x[2]
        n_threads = x[3]
        """

        print(x)
        tpr = x[0]
        out = x[1]

        # run simulation
        cmd = "gmx mdrun -v -s %s -deffnm %s"%(tpr, out)
        print(" running: ",cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        return (out, err)

    def run_md(self, n_threads, params):
        """
        worker provisional para 'Equilibration'
        """

        # note for the moment wet_xtc_file, dry_xtc_file and last_wet_snapshot
        # overwrite each other, until I see their need

        cmd = 'gmx grompp -f %s -c processed_%s -p topol.top -maxwarn 1'\
                % (params.filemdp, self.pdb); \
                'gmx mdrun -nt %d -s topol.tpr -x %s -c processed_%s -g prod.log' % \
                (n_threads, self.pdb, self.wet_xtc_file)
        print(" running: ",cmd)
        p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()

        # if running in production mode, trjconv to create a dry XTC
        if params.md_step == 'Production':
            assert self.dry_xtc_file is not None and self.last_wet_snapshot is not None
            cmd = 'echo 0 | gmx trjconv -f end.gro -s topol.tpr -o %s -pbc whole' % self.last_wet_snapshot
            print(" running: ",cmd)
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            #print(output, error)
            cmd = 'echo PROTEIN | gmx trjconv -f %s -s topol.tpr -o %s -pbc whole' % (self.wet_xtc_file, self.dry_xtc_file)
            print(" running: ",cmd)
            p = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = p.communicate()
            #print(output, error)

        return (out, err)

    def clean_working_directory(self):

        print ("Cleaning up...")
        os.system('rm \#*')
        os.system('mkdir logs')
        os.system('mv *.log logs')
        os.system('rm *.trr *.top *.itp *.edr *.cpt *.tpr out.gro conf.gro')
        os.system('mv equilibration.xtc *.mdp *.gro logs')
    
        return