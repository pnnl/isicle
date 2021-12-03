import isicle
import os
import sphinx_skeleton as ss

if __name__ == "__main__":
    # generic package reference
    package = isicle
    
    # modules
    modules = ss.get_modules(os.path.dirname(package.__file__))

    # documentation layout
    doc_layout = {'getting_started': ['installation',
                                      'tutorial'],
                  # edit user_guide per package
                  'user_guide': ['loading_saving',
                                 'molecule_prep',
                                 'conformer_generation',
                                 'structure_optimization',
                                 'property_prediction'],
                  'api_reference': [package.__name__] + modules,
                  'project_info': ['faq',
                                   'citing_and_citations',
                                   'contributing',
                                   'license',
                                   'disclaimer']}

    # readable titles
    titles = {'api_reference': 'API reference',
              'getting_started': 'Getting started',
              'project_info': 'Project info',
              'user_guide': 'User guide',
              'installation': 'Installation',
              'tutorial': 'Tutorial',
              'loading_saving': 'Loading/saving',
              'molecule_prep': 'Molecule preparation',
              'conformer_generation': 'Conformer generation',
              'structure_optimization': 'Structure optimization',
              'property_prediction': 'Property prediction',
              'citing_and_citations': 'Citing and citations',
              'faq': 'Frequently asked questions',
              'contributing': 'Contributing',
              'license': 'License',
              'disclaimer': 'Disclaimer'}

    # license
    with open(os.path.join(os.path.dirname(os.path.dirname(package.__file__)), 'LICENSE')) as f:
        license = f.read()

    # instance
    skel = ss.SphinxSkeletonizer(package)
    skel.set_tree(doc_layout)
    skel.set_titles(titles)
    skel.clean()

    # custom write functions
    def api_reference_fxn(key):
        if package.__name__ not in key:
            fullkey = '{}.{}'.format(package.__name__, key)
        else:
            fullkey = key

        string = ss.generate_title(key)
        string += '.. automodule:: {}\n'.format(fullkey)
        string += '\t:members:\n'
        string += '\t:private-members:\n'
        string += '\t:undoc-members:\n'

        # top level
        if package.__name__ in key:
            string += '\t:imported-members:\n'
        
        return string

    
    def license_fxn(key):
        string = ss.generate_title(key)
        string += license
        return string

    
    def citation_fxn(key):
        string = ss.generate_title(key)
        string += ('If you would like to reference {} in an academic paper, '
                   'we ask you include the following:\n'
                   '[citations here]\n').format(package.__name__)
        return string

    # set write functions
    for page in doc_layout['api_reference']:
        skel.set_template(page, api_reference_fxn)

    skel.set_template('license', license_fxn)
    skel.set_template('citing_and_citations', citation_fxn)

    skel.build_tree()
    skel.write_index()
