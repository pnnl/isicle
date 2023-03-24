from openbabel import pybel
from openbabel import openbabel
import numpy as np
import logging


# Class for representing box
class Box:
    """
    This class represents the dimensions based on the molecule geometry
    """

    def __init__(self, xyzr, h=0.5, ext=4):
        self.logger = logging.getLogger(__name__)
        self.xmin = []
        self.xmax = []
        self.ymin = []
        self.ymax = []
        self.zmin = []
        self.zmax = []
        self.h = h
        self.ext = ext
        self.nx = []
        self.ny = []
        self.nz = []
        self.x = []
        self.y = []
        self.z = []
        self.xyzr = xyzr

        """
        This function sets the box dimensions based on the molecule geometry.
        INPUT:
            h: grid spacing in x,y, and z directions of the mesh on to
                         which the molecule geometry is mapped.
            xyzr: array of x,y,z, and radius values for each atom in the molecule
            ext:  distance of the box edge from the atom closest to the edge.

        OUTPUT:
            pybelmol: the new pybel molecule object
            total_chg: total charge on the new molecule
            mol3Dfile:  the .mol file of the adduct
        """
        # Determine the edge positions of the box along the x-axis
        [xleft, xright] = self.getBoxEdges(0)

        # Determine the edge positions of the box along the x-axis
        [yleft, yright] = self.getBoxEdges(1)

        # Determine the edge positions of the box along the x-axis
        [zleft, zright] = self.getBoxEdges(2)

        nx = np.int(np.round((xright - xleft) / self.h) + 1)
        ny = np.int(np.round((yright - yleft) / self.h) + 1)
        nz = np.int(np.round((zright - zleft) / self.h) + 1)

        self.logger.debug("nx = %s ny = %s nz = %s", nx, ny, nz)

        # Correct the values of edge distance based on nx, ny, and nz values.
        xright = xleft + (nx - 1) * self.h
        yright = yleft + (ny - 1) * self.h
        zright = zleft + (nz - 1) * self.h

        self.xleft = xleft
        self.xright = xright
        self.yleft = yleft
        self.yright = yright
        self.zleft = zleft
        self.zright = zright

        self.nx = nx
        self.ny = ny
        self.nz = nz

        for i in range(0, nx):
            xi = self.xleft + i * self.h
            self.x.append(xi)

        for i in range(0, ny):
            yi = self.yleft + i * self.h
            self.y.append(yi)

        for i in range(0, nz):
            zi = self.zleft + i * self.h
            self.z.append(zi)

    # Function to determine the box length dimension along an axis direction
    def getBoxEdges(self, index):
        # Determine the edge positions of the box along the axis
        natoms = len(self.xyzr)
        patom_min = self.xyzr[0][index] - self.xyzr[0][3] - self.ext
        patom_max = self.xyzr[0][index] + self.xyzr[0][3] + self.ext

        for i in range(1, natoms):  # i:1, 2, ...,natoms - 1
            tmp = (
                self.xyzr[i][index] + self.xyzr[i][3] + self.ext
            )  # tmp = p + rad + ext
            if patom_max < tmp:
                patom_max = tmp
            tmp = (
                self.xyzr[i][index] - self.xyzr[i][3] - self.ext
            )  # tmp = p - rad - ext
            if patom_min > tmp:
                patom_min = tmp

        left_ex = patom_min / self.h
        right_ex = patom_max / self.h

        pleft = np.floor(left_ex) * self.h - self.ext
        pright = np.ceil(right_ex) * self.h + self.ext
        return pleft, pright

    # create atom accessibility map
    def createAtomAccessibiltyMap(self, prob_rad):
        """
        This function creates the molecule surface accessibility map of the new
        atom to be added in the molecule.
        INPUT:
            prob_rad: radius of the atom to be added

        OUTPUT:
            pybelmol: the new pybel molecule object
            total_chg: total charge on the new molecule
            mol3Dfile:  the .mol file of the adduct
        """

        # initialize the values of the accessibilty map to zero
        acc_map = np.zeros((self.nx, self.ny, self.nz))

        natoms = len(self.xyzr)

        # Create the  accessibility map by
        # rolling a spherical particle of radius (= probe radius) on the surface of each atom
        # to determine the domain accessible/in-accessible to the particle.
        for i in range(0, natoms):
            atom_radius = self.xyzr[i][3]
            atom_x = self.xyzr[i][0]
            atom_y = self.xyzr[i][1]
            atom_z = self.xyzr[i][2]

            r_atom = atom_radius + prob_rad

            x1 = np.floor((atom_x - self.xleft - r_atom) / self.h)
            x2 = np.ceil((atom_x - self.xleft + r_atom) / self.h)
            ix1 = np.int(x1)
            ix2 = np.int(x2)

            y1 = np.floor((atom_y - self.yleft - r_atom) / self.h)
            y2 = np.ceil((atom_y - self.yleft + r_atom) / self.h)
            iy1 = np.int(y1)
            iy2 = np.int(y2)

            z1 = np.floor((atom_z - self.zleft - r_atom) / self.h)
            z2 = np.ceil((atom_z - self.zleft + r_atom) / self.h)
            iz1 = np.int(z1)
            iz2 = np.int(z2)

            for ix in range(ix1, ix2 + 1):
                for iy in range(iy1, iy2 + 1):
                    for iz in range(iz1, iz2 + 1):
                        xdist = self.x[ix] - atom_x
                        ydist = self.y[iy] - atom_y
                        zdist = self.z[iz] - atom_z
                        rdist = np.sqrt(xdist**2 + ydist**2 + zdist**2)

                        if rdist < r_atom:
                            acc_map[ix, iy, iz] = 1.0
        return acc_map


def addAtomToMol(mol, atom, idx, covalent=True):
    """
    This function adds an atom (atom) to a molecule (mol) at the site
    near an atom (idx). It uses the covalent radius of the atoms
    to determine the accessible site for placing the the new atom.

    INPUT:
        mol - pybel molecule object
        atom - element symbol (string) of the atom to be added
        idx - atom number of the atom near which the new atom is
              placed
        covalent - True for covalent; otherwise non-covalent (van der Waal)
    OUTPUT:
        mol - pybel molecule object containing the new atom
    """

    logger = logging.getLogger(__name__)
    # openbabel 2.3.1 approach
    # etab = openbabel.OBElementTable()  # element object
    # openbabel 3.0.1 approach
    # etab functions moved to openbabel.openbabel namespace

    atoms_radius = []
    xyzr = []

    for x in mol.atoms:
        atomic_num = x.atomicnum  # get atomic number of each atom
        if covalent is True:
            # use covalent radius
            atom_rad = openbabel.GetCovalentRad(atomic_num)

        else:
            # use van der Waals radius (but for now use covalent)
            atom_rad = openbabel.GetVdwRad(atomic_num)
        atoms_radius.append(atom_rad)

        a = x.coords + (atom_rad,)  # create x,y,z,radius array
        xyzr.append(a)

    # get the atomic number and covalent radius of the new atom to be added
    newatom_atomic_num = openbabel.__dict__[atom]

    if covalent is True:
        newatom_rad = openbabel.GetVdwRad(newatom_atomic_num)
    else:
        newatom_rad = openbabel.GetVdwRad(1)
    newatom_xyzr = [0, 0, 0, newatom_rad]

    # set box dimensions based on molecule size
    box = Box(xyzr, 0.5, 10)

    # create accessibilty map
    probe_radius = newatom_rad
    acc_map = box.createAtomAccessibiltyMap(probe_radius)
    logger.info("Finished creating atom accessibility map.")
    xc = xyzr[idx][0]
    yc = xyzr[idx][1]
    zc = xyzr[idx][2]
    rc = xyzr[idx][3]

    rdist = rc + probe_radius
    logger.info("rdist = %s", rdist)

    # randomly select a position for the new atom near the specified atom of the
    # molecule
    found = 0
    delr = 0.01
    ntrials = 5000
    ndist = 50
    count = 0

    logger.info(
        "Selecting a position for the new atom near the atom located at: (%s, %s, %s)",
        xc,
        yc,
        zc,
    )

    while found == 0 and count < ndist:
        logger.debug("count #%s", count)
        rdist = rdist + delr
        ntrials = 5000
        trial = 0

        x1 = xc - rdist
        x2 = xc + rdist
        if x1 < box.xleft:
            x1 = box.xleft

        if x2 > box.xright:
            x2 = box.xright

        y1 = yc - rdist
        y2 = yc + rdist
        if y1 < box.yleft:
            y1 = box.yleft

        if y2 > box.yright:
            y2 = box.yright

        z1 = zc - rdist
        z2 = zc + rdist
        if z1 < box.zleft:
            z1 = box.zleft

        if z2 > box.zright:
            z2 = box.zright

        logger.debug("x1, y1, z1: %s %s %s", x1, y1, z1)
        logger.debug("x2, y2, z2: %s %s %s", x2, y2, z2)
        while trial < ntrials and found == 0:
            logger.debug("trial #: %s", trial)
            rand_x = x1 + np.random.ranf() * (x2 - x1)
            rand_y = y1 + np.random.ranf() * (y2 - y1)
            rand_z = z1 + np.random.ranf() * (z2 - z1)

            ix = np.int((rand_x - box.xleft) / box.h)
            iy = np.int((rand_y - box.yleft) / box.h)
            iz = np.int((rand_z - box.zleft) / box.h)

            if acc_map[ix, iy, iz] == 0:
                newatom_xyzr[0] = rand_x
                newatom_xyzr[1] = rand_y
                newatom_xyzr[2] = rand_z
                found = 1
                logger.info("Found a site: (%s, %s, %s)", rand_x, rand_y, rand_z)
                rd = np.sqrt(
                    (rand_x - xc) * (rand_x - xc)
                    + (rand_y - yc) * (rand_y - yc)
                    + (rand_z - zc) * (rand_z - zc)
                )
                logger.info("inter-atom distance (A): %s", rd)

            else:
                trial = trial + 1
                found = 0

        count = count + 1

    if found == 0:
        logger.critical("Could not find a site to place adduct")
        raise Exception("Could not find a site to place adduct")
    else:
        newatom = mol.OBMol.NewAtom()
        newatom.SetAtomicNum(newatom_atomic_num)
        newatom.SetVector(rand_x, rand_y, rand_z)

        xyzr.append(newatom_xyzr)

        if covalent is True:
            new_atomnum = len(mol.atoms)
            mol.OBMol.AddBond(idx, new_atomnum, 1)
            # mol.OBMol.AddBond(idx + 1, new_atomnum, 1)

            # good for single bonds. For double and triple bonds, we have to
            # the appropriate covalent radius. For our purposes, we only need
            # single bond (for protonation).

    return mol


def removeAtomFromMol(mol, idx):
    """
    This function removes an atom from the pybel moleule object.
    """

    # mol should contain the atom coordinates
    # atom = mol.OBMol.GetAtomById(idx)
    atom = mol.OBMol.GetAtom(idx)
    # delete the atom
    mol.OBMol.DeleteAtom(atom)

    return mol


def nearestHydrogen(mol, idx):
    logger = logging.getLogger(__name__)
    # iatom = mol.atoms[idx]
    iatom = mol.OBMol.GetAtom(idx)
    logger.debug("Starting atom: %s, type %s", idx, iatom.GetAtomicNum())

    # get the neighboring atoms of the selected atom
    nbatoms = openbabel.OBAtomAtomIter(iatom)
    for nb in nbatoms:
        logger.debug("Connected to: %s, type %s", nb.GetId(), nb.GetAtomicNum())
        # if the neighboring atom is hydrogen
        if nb.GetAtomicNum() == 1:
            return nb.GetIdx()

    raise Exception("Hydrogen not found.")
