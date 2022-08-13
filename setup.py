from distutils.core import setup


install_requires = []
with open("requirements.txt") as f:
    install_requires = f.read().splitlines()


setup(
    name = 'SolOrder',
    packages = ['SolOrder'],
    version = '1.0.9',  # Ideally should be same as your GitHub release tag varsion
    description = 'Updated a few..',
    author = 'Dibyendu Maity',
    author_email = 'dibyendumaity1999@bose.res.in',
    url = 'https://github.com/dmighty007/SolOrder',
    download_url = 'https://github.com/dmighty007/SolOrder',
    keywords = ['water', 'SOL', 'Order','Parameter'],
    classifiers = [],
    project_urls = {
        "Bug Tracker": "https://github.com/dmighty007/SolOrder/issues"
    },
    license='MIT',
    install_requires=install_requires,
)
