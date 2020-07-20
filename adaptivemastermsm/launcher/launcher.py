"""
This file is part of the AdaptiveMasterMSM package.

"""
# General
import os
import subprocess
from cmds import cmds

class launcher(object):
    """
    Call GROMACS and process output to pass it to class analyzer
    """

    def __init__(self, sysfile, md_step, n_threads, forcefield, \
            pdb_fn, wet_xtc_file,  dry_xtc_fn=None, last_wet_snapshot_fn=None):
        """
        Read input from system, run MD, and process output
        for the class analyzer

        Args:
            sysfile (str): File containing system information
            md_step (str): 'Equilibration' or 'Production'
            forcefield (str): amber96, ...
            pdb_fn (str): PDB filename
            wet_xtc_file (str)
            dry_xtc_file (str)
            last_wet_snapshot_file (str)

        Return???:
            trajfile (str): File containing processed trajectories
        """

        print "Running on %d threads" % n_threads
        # create an instance of system class
        p = system(water, md_step)

        # build a list of commands to execute, then shell them out
        cmds = []

        if p.md_step == 'Production':
            cmds.append('pdb2gmx -f %s -ff %s -water %s' % (pdb_fn, forcefield, p.water))
            # feed pdb2gmx out straight into grompp
            cmds.append('cp conf.gro start.gro')
            cmds.extend( run_md(n_threads, p) )

        elif p.md_step == 'Equilibration':
            # first equilibrate
            cmds.append('pdb2gmx -f %s -ff %s -water %s -ignh' % (pdb_fn, forcefield, p.water))
            cmds.extend( p.build_box() )
            user_wet_xtc_file = wet_xtc_file
            self.wet_xtc_file = 'equilibration.xtc'
            cmds.extend( run_md(n_threads, p) )

            # then production
            #p.md_step = 'Production'  #ionix, ojo aqui
            #p.init_parameters() actualiza self.run y crea otro md_run
            p_prod = system(water, 'Production')
            cmds.append('cp end.gro start.gro')
            self.wet_xtc_file = user_wet_xtc_file
            cmds.extend( run_md(n_threads, p_prod) )

        shell_out(cmds)
        clean_working_directory()

        return

    def run_md(n_threads, Parameters):
    
        #write_mdp('%s_parameters.mdp' % Parameters.md_step)
        cmds = ['grompp -f %s_parameters.mdp -c start.gro -p topol.top -maxwarn 1'\
                % Parameters.md_step, 'mdrun -nt %d -s topol.tpr -x %s -c end.gro -g prod.log'\
             % (self.n_threads, self.wet_xtc_fn) ]

        # if running in production mode, trjconv to create a dry XTC
        if Parameters.md_step == 'Production':
            assert self.dry_xtc_fn is not None and self.last_wet_snapshot_fn is not None
            cmds.append('echo 0 | trjconv -f end.gro -s topol.tpr -o %s -pbc whole' % self.last_wet_snapshot_fn)
            cmds.append('echo PROTEIN | trjconv -f %s -s topol.tpr -o %s -pbc whole' % (self.wet_xtc_fn, self.dry_xtc_fn))

        return cmds


    def shell_out(commands):
    """ executes commands intelligently via subprocess -
         waits for serial completion
    """
    
        log = open('driver.log','w')
        for command in commands:
            print "Executing: %s" % command
            x = subprocess.Popen( command, bufsize=-1, shell=True, stdout=log, stderr=log )
            x.wait()
            if x.returncode != 0:
                raise Exception("Error: '%s' quit with exit status: %d" % (command, x.returncode) )
        log.close()
    
        return

    def clean_working_directory():

        print "Cleaning up..."
        os.system('rm \#*')
        os.system('mkdir logs')
        os.system('mv *.log logs')
        os.system('rm *.trr *.top *.itp *.edr *.cpt *.tpr out.gro conf.gro')
        os.system('mv equilibration.xtc *.mdp *.gro logs')
    
        return