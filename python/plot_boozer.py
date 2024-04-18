import matplotlib.pyplot as plt
from libneo.boozer import BoozerFile

booz = BoozerFile("in_file")
contours = booz.get_contours_in_r_z_plane(0.0, 10)

plt.figure(figsize=(10, 10))
plt.plot(contours[0], contours[1], 'k,')
