"""
This file is part of the AdaptiveMasterMSM package.

"""

class System(object):
    """
    Create files to be able to run GROMACS
    """

    def __init__(self, water, md_step):
        """
        Args:
            water: water model for GROMACS (tip3p, ...)
            md_step (str): must be 'Equilibration' or 'Production'

        """
        # Read user-defined parameters
        self.water = water
        self.md_step = md_step

        # Take production or equilibration paths
        if self.md_step == 'Equilibration':
            self.driver_equilibration()
        elif self.md_step == 'Production':
            self.driver_production()
        else:
            print ("md_step %s not valid" % self.md_step)
            raise Exception("It must be 'Production' or 'Equilibration'")

        self.write_mdp('%s_parameters.mdp' % self.md_step)

        return
        
    def driver_equilibration(self):
    
        self.equil = self.set_equilibration()
        self.run = {}
        self.equil['total_steps'] = \
            int(self.equil['total_time'] / self.equil['timestep'])
        
        for p in self.equil.keys():
            self.run[p] = self.equil[p]

        return           

    def driver_production(self):

        self.prod = self.set_production()
        self.run = {}
        self.prod['total_steps'] = \
            int(self.prod['total_time'] / self.prod['timestep'])

        for p in self.prod.keys():
            self.run[p] = self.prod[p]

        return

    def set_production(self):
        """
        Parameters for the MD production run
        """
        production = {
            'timestep'       : 0.002,         # ps
            'total_time'     : 1,             # ps / 10ns
            'log_freq'       : 10000,         # timesteps / 20ps
            'xtc_freq'       : 100,         # timesteps
            'temperature'    : 300,           # Kelvin
            'thermostat'     : 'Nose-Hoover', # type
            'box_size'       : 3.5,           # Angstroms
            'barostat'       : 'No',          # No barostat
            'pressure'       : 1,             # pressure
            'Cl'             : 0,             # number Cl's
            'Na'             : 0              # number Na's
        }

        return production

    def set_equilibration(self):
        """
        Parameters for the MD equilibration run
        """
        equilibration = {
            'timestep'       : 0.002,         # ps
            'total_time'     : 1,          # ps
            'barostat'       : 'berendsen',   # type
            'pressure'       : 1,              # bar
            'log_freq'       : 10000,         # timesteps / 20ps
            'xtc_freq'       : 100,         # timesteps
            'temperature'    : 300,           # Kelvin
            'thermostat'     : 'Nose-Hoover', # type
            'box_size'       : 3.5,           # Angstroms
            'Cl'             : 0,             # number Cl's
            'Na'             : 0              # number Na's
        }

        return equilibration

    def build_box(self):

        self.write_minimization_mdp('Minimization.mdp')

        # choose a water topology file
        water_dict  = {'tip3p' : 'spc216.gro'}
        if self.water in water_dict.keys():
            water_topol = water_dict[self.water]
        else:
            print ("Could not find a topology in 'w_top_dict' for water: %s" % self.water)
            raise Exception("Invalid water model for GROMACS.")

        # format the water string
        ion_str = ''
        if self.run['Cl'] > 0:
            ion_str += '-nn %d ' % self.run['Cl']
        if self.run['Na'] > 0:
            ion_str += '-np %d ' % self.run['Na']

        cmds =[
        'gmx editconf -f conf.gro -bt cubic -box %f %f %f -align 1 1 1' %((self.run['box_size'],)*3),
        'gmx genbox -cp out.gro -cs %s -p topol.top' % water_topol,
        'gmx grompp -f Minimization.mdp -c out.gro -p topol.top',
        'echo SOL | gmx genion  -s topol.tpr -o out.gro -p topol.top %s' % ion_str,
        'gmx grompp -f Minimization.mdp -c out.gro -p topol.top',
        'gmx mdrun -s topol.tpr -x minimization.xtc -c start.gro -pd -g EM.log'
        ]

        return cmds
   
    def write_mdp(self, f_name):

        txt = """; GENERATED BY MONAKOS: gromacs_driver.py
    integrator               = md
    dt                       = %f
    nsteps                   = %d
    nstxout                  = 0
    nstvout                  = 0
    nstlog                   = %d
    nstenergy                = 0
    nstxtcout                = %d
    xtc_grps                 = System
    nstlist                  = 0
    ns_type                  = grid
    pbc                      = xyz
    periodic_molecules       = no
    rlist                    = 0.03
    coulombtype              = PME
    fourier_nx               = 0
    fourier_ny               = 0
    fourier_nz               = 0
    optimize_fft             = yes
    pme_order                = 4
    fourierspacing           = 0.08
    rcoulomb                 = 0.1
    vdwtype                  = shift
    rvdw                     = 0.1
    ; rvdw_switch              = 1.0
    tcoupl                   = %s
    tc_grps                  = System
    tau_t                    = 1
    ref_t                    = %f
    Pcoupl                   = %s
    tau_p                    = 1
    pcoupltype               = isotropic
    compressibility          = 0.000045
    ref_p                    = %f
    gen_vel                  = yes
    gen_temp                 = %f
    constraints              = hbonds
    continuation             = no
    ; continuation             = yes
    morse                    = no
    implicit_solvent         = no
    """ % ( self.run['timestep'],
        self.run['total_steps'],
        self.run['log_freq'],
        self.run['xtc_freq'],
        self.run['thermostat'],
        self.run['temperature'],
        self.run['barostat'],
        self.run['pressure'],
        self.run['temperature'] )

        f = open(f_name, 'w')
        f.write(txt)
        f.close()
        print ("Generated: %s" % f_name)

        return

    def write_minimization_mdp(self, f_name):

        txt = """; GENERATED BY MONAKOS: gromacs_driver.py
    integrator               = steep
    emstep                   = 0.001
    emtol                    = 10.0
    nsteps                   = 100000
    nstxout                  = 0
    nstvout                  = 0
    nstlog                   = 10
    nstenergy                = 0
    nstxtcout                = 0
    xtc_grps                 = System
    energygrps               = 
    nstlist                  = 1
    ns_type                  = grid
    pbc                      = xyz
    periodic_molecules       = no
    rlist                    = 1.5
    rcoulomb                 = 1.5
    rvdw                     = 1.5
    tcoupl                   = no
    Pcoupl                   = no
    gen_vel                  = no
    constraints              = none
    continuation             = yes
    morse                    = no
    implicit_solvent         = no
    """
        f = open(f_name, 'w')
        f.write(txt)
        f.close()
        print ("Generated: %s" % f_name)

        return


