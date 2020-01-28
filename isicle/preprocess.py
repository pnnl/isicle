import subprocess
import multiprocessing as mp


def inchi2smi(inchi):
    '''Converts InChI string to SMILES string.'''

    try:
        res = subprocess.check_output('obabel -iinchi -:"%s" -ocan' % inchi,
                                      stderr=subprocess.DEVNULL, shell=True).decode('ascii')
    except:
        print(inchi, 'failed.')
        return None

    return res.strip()


def smi2inchi(smi):
    '''Converts SMILES string to InChI string.'''

    try:
        res = subprocess.check_output('obabel -ismi -:"%s" -oinchi' % smi,
                                      stderr=subprocess.DEVNULL, shell=True).decode('ascii')
    except:
        print(smi, 'failed.')
        return None

    return res.strip()


def canonicalize(smi):
    '''Canonicalizes SMILES string.'''

    try:
        res = subprocess.check_output('obabel -ismi -:"%s" -ocan' % smi,
                                      stderr=subprocess.DEVNULL, shell=True).decode('ascii')
    except:
        print(smi, 'failed.')
        return None

    return res.strip()


def desalt(smi):
    try:
        res = subprocess.check_output('obabel -ismi -:"%s" -r -ocan' % smi,
                                      stderr=subprocess.DEVNULL, shell=True).decode('ascii')
    except:
        print(smi, 'failed.')
        return None

    return res.strip()


def neutralize(smi):
    '''Neutralizes an canonical SMILES string (alternate).'''

    def neutralize_inchi(inchi):
        '''Neutralizes an InChI string.'''
        if 'q' in inchi:
            layers = inchi.split('/')
            new = layers[0]
            for i in range(1, len(layers)):
                if 'q' not in layers[i]:
                    new += '/' + layers[i]
            return new
        return inchi

    inchi = smi2inchi(smi)
    inchi = neutralize_inchi(inchi)
    return inchi2smi(inchi)


def tautomerize(smi):
    try:
        res = subprocess.check_output('cxcalc majortautomer -f smiles %s' % smi,
                                      stderr=subprocess.DEVNULL, shell=True).decode('ascii')
    except:
        print(smi, 'failed.')
        return None

    return res.strip()


def _process(smiles):
    return canonicalize(tautomerize(neutralize(desalt(smiles))))


def process(smiles):
    with mp.Pool(processes=mp.cpu_count()) as p:
        return p.map(_process, smiles)
