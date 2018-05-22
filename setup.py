from setuptools import setup, find_packages, findall

setup(  name                = 'scud',
        version             = '0.0.0',
        description         = 'diffuse scattering tools',
        author              = 'T de Klijn',
        author_email        = 'timdeklijn@gmail.com',
        url                 = '',
        license             = '',
        install_requires    = ( ),
        package_dir         = {'':'scud'},
#        scripts             = install_scripts,
        packages            = find_packages(where   ='scudlib-python',
                                            exclude = ( )
                                            ),
        include_package_data= True,
        classifiers         = (  ),
        )
