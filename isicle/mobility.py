import isicle
from isicle.interfaces import WrapperInterface
import os
import subprocess
import tempfile
from pkg_resources import resource_filename
import shutil


def calculate_ccs(geom, **kwargs):
    # Initialize wrapper
    mw = MobcalWrapper()

    # Set geometry
    mw.set_geometry(geom)

    # Save geometry
    mw.save_geometry()

    # Configure
    mw.configure(**kwargs)

    # Save config
    mw.save_config()  # Currently does nothing...

    # Run mobility calculation
    mw.run()

    # Finish/clean up
    res = mw.finish()

    # TODO: update geometry object with CCS result

    return res


class MobcalWrapper(WrapperInterface):

    def __init__(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.atom_params = os.path.join(self.temp_dir.name,
                                        'atomtype_parameters.in')
        self.mobcal_params = os.path.join(self.temp_dir.name,
                                          'mobcal.params')

        self.infile = os.path.join(self.temp_dir.name,
                                   self.geom.basename + '.mfj')
        self.outfile = os.path.join(self.temp_dir.name, self.geom.basename + '.out')
        self.logfile = os.path.join(self.temp_dir.name, self.geom.basename + '.log')

    def set_geometry(self, geom):
        '''
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        '''

        # Assign geometry
        self.geom = geom

    def save_geometry(self):
        '''
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        Raises
        ------
        TypeError
            If geometry loaded from .xyz is saved to another format.

        '''

        # All other formats
        self.geom.save(self.infile)

    def _configure_lennard_jones(self, path=None):
        if path is None:
            path = resource_filename('isicle', 'resources/lennard_jones.txt')

        shutil.copy2(path, self.atom_params)

    def _configure_mobcal(self, i2=5013489, buffer_gas='helium',
                          buffer_gas_mass=4.0026, temp=300, ipr=1000,
                          itn=10, inp=48, imp=1024, processes=24):

        d = {'I2': i2,
             'BUFFER_GAS': buffer_gas.upper(),
             'BUFFER_GAS_MASS': buffer_gas_mass,
             'TEMPERATURE': temp,
             'IPR': ipr,
             'ITN': itn,
             'INP': inp,
             'IMP': imp,
             'NUM_THREADS': processes}

        with open(self.mobcal_params, 'w') as f:
            f.write('\n'.join(['{} {}'.format(k, v) for k, v in d.items()]))

    def configure(self, lennard_jones='default', i2=5013489,
                  buffer_gas='helium', buffer_gas_mass=4.0026, temp=300,
                  ipr=1000, itn=10, inp=48, imp=1024, processes=24):

        # Handle default case
        if lennard_jones == 'default':
            lennard_jones = None

        # Configure Lennard-Jones potentials
        self._configure_lennard_jones(lennard_jones)

        # Configure Mobcal parameters
        self._configure_mobcal(i2=i2, buffer_gas=buffer_gas,
                               buffer_gas_mass=buffer_gas_mass,
                               temp=temp, ipr=ipr, itn=itn, inp=inp, imp=imp,
                               processes=processes)

    def save_config(self):
        pass

    def run(self):
        subprocess.call('mobcal {} {} {} {} &> {}'.format(self.mobcal_params,
                                                          self.atom_params,
                                                          self.infile,
                                                          self.outfile,
                                                          self.logfile),
                        shell=True)

    def finish(self):
        # Initialize parser
        parser = isicle.parse.MobcalParser()

        # Load output file
        parser.load(self.outfile)

        # Extract result
        result = parser.parse()

        # Remove temporary files
        self.temp_dir.cleanup()

        return result
