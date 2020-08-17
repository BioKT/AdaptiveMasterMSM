from adaptivemastermsm.launcher import launcher
###from adaptivemastermsm.system import system
lau = launcher.Launcher('Equilibration', 4, 'amber96', 'tip3p', 'alaTB.gro', 'wetfn', dry_xtc_file='dryfn', last_wet_snapshot='lastwetfn')

