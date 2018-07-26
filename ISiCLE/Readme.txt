#
#
#
#######################################################################
#######################################################################
###																	###
###                             ISICLE							    ###
###																	###
#######################################################################
#######################################################################
#
#
#
# Type the command at the bottom of this file in a MS-DOS command line 
# in Microsoft Windows 
# 
#
#
# Replace every <network-id> with your Cascade user account id
#
#
# Replace every <account-id> with your Cascade project account id 
#
#
# Replace <Molecule file with the Microsoft Excel extension> with your 
# Excel file name
# Ex: MoleculeSet.xlsx
#
#
# Replace <Method file with the Microsoft Excel extension> with your 
# Excel file name
# Ex: Method.xlsx
#
#
# Replace <nnodes> with the number of nodes you will use
# Ex: 6
# 
#
# Replace <nprocessors> with the number of processors you will use
# Ex: 16
#
#
# Replace <walltime> with the wall time you will use
# Ex: 00:15:00
#
#

python isicle.py -accnt <account-id> -molcs <Molecule file with Microsoft Excel extension> -metd <Method file with Microsoft Excel extension> -host <network-id>@cascade.emsl.pnl.gov -dest /dtemp/<network-id>/ -nn <nnodes> -nc <nprocessors> -time <walltime>

Ex: python isicle.py -accnt <account-id> -molcs MoleculeSet.xlsx -metd Method.xlsx -host <network-id>@cascade.emsl.pnl.gov -dest /dtemp/<account-id>/<network-id>/ -nn 6 -nc 16 -time 00:15:00
