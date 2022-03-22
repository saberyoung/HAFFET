import json
import numpy as np



class CatGen(object):

    def __init__(self, transient_name):
        self.transient_name = transient_name
        self.photometry = []
        self.lumdist = None
        self.lumdist_err = None

    def dump_json(self, destfile):
        with open(destfile, 'w') as f:
            ast = self.dump_ast()
            #for x in self.photometry:
            #    print(repr(x))
            #    json.dumps(x)
            json.dump(ast, f)
        
    def dump_ast(self):
        evt = {
            'name': self.transient_name,
            'sources': [
                {'name': 'gp3', 'alias': '1'}
            ],
            'alias': [
                {'value': self.transient_name, 'source': '1'}
            ],
            'photometry': self.photometry
        }

        if self.lumdist is not None:
            evt['lumdist'] = [
                {
                    "value": '%e' % self.lumdist,
                    "u_value": "Mpc",
                    "source": "1",
                }
            ]
            if self.lumdist_err is not None:
                evt['lumdist'][0]['e_value'] = '%e' % self.lumdist_err
        # TODO: check for redshift, etc. and add that..

        ret = {
            self.transient_name: evt
        }
        return ret

    def add_photometry_point(self, mjd, telescope, instrument, band, mag, mag_err, mag_system):
        p = {
            'time': float(mjd),
            'band': band,
            'instrument': instrument,
            'telescope': telescope,
            'magnitude': float(mag),
            'e_magnitude': float(mag_err),
            'u_time': 'MJD',
            'source': '1',
            'system': mag_system
        }
        if telescope is None:
            del p['telescope']
        if instrument is None:
            del p['instrument']
        if mag_err is None:
            del p['e_magnitude']
            p['upperlimit'] = True
        self.photometry.append(p)

    def add_lumdist(self, lumdist, lumdist_err=None):
        assert self.lumdist is None
        self.lumdist = lumdist
        self.lumdist_err = lumdist_err
        


def convert_mosfit(df, target_name, mag_column='MAG', magerr_column='MAGERR', upper_lim_mag_column='LIMIT', mag_system='AB'):
    # construct photometry table
    photometry = []

    # get rid of multi-index and convert them to columns
    df = df.reset_index()
    for idx, row in df.iterrows():
        p = {
            'time': row['MJD'],
            'band': row['filter'],
            'instrument': row['instrument'],
            'telescope': row['telescope'],
            'magnitude': row[mag_column],
            'e_magnitude': row[magerr_column],
            'u_time': 'MJD',
            'source': '1',
            'system': mag_system
        }
        if p['telescope'] is None:
            del p['telescope']
        if p['instrument'] is None:
            del p['instrument']
        if np.isnan(row[mag_column]):
            if upper_lim_mag_column is None:
                continue
            
            p['magnitude'] = row[upper_lim_mag_column]
            del p['e_magnitude']
            p['upperlimit'] = True
            
        
        photometry.append(p)
        

    evt = {
        'name': target_name,
        'sources': [
            {'name': 'gp3', 'alias': '1'}
        ],
        'alias': [
            {'name': target_name, 'source': '1'}
        ],
        'photometry': photometry
    }

    ret = {
        target_name: evt
    }
    return ret

