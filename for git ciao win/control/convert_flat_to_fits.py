import numpy as np
import astropy.io.fits as fits

flatfile = open(r'C:\Users\AOC\bin\zygo_flat.txt', 'r')
vLines = flatfile.readlines()

flat = np.zeros(len(vLines)).astype(np.float32)

count = 0
for line in vLines:
    value = np.float32(line.rstrip('\n'))
    print('value >> ', value)
    flat[count] = value
    count += 1

filename = r'C:\Users\AOC\bin\zygo_flat.fits'
fits.writeto(filename,flat,overwrite=True)
fits.setval(filename,'ISFLAT',value=1)

print('Done.')
# print(vLines)
# print(len(vLines))
