This is the README for the model associated with the paper:

Birdno MJ, Kuncel AM, Dorval AD, Turner DA, Gross RE, Grill WM (2012)
Stimulus features underlying reduced tremor suppression with
temporally patterned deep brain stimulation. J Neurophysiol 107:364-83

Thalamic Network Model
Merrill J. Birdno
May 2009

TC Neuron Derived from McIntyre et.el., 2004

BASIC INSTRUCTIONS:

1) You will need to use mknrn.dll (see NEURON documentation)(MAC or
mswin) or nrnivmodl (unix/linux) to create a new 'special' file with
the compiled .mod mechanisms. Copy this special file into the 'master'
simulation directory included with the model.

2) ~Line 627 in DBSstim.hoc, you will see the line:

		indrt = strobj.tail(subdir,"ModelDirectory/","x")
	
	You will need to change "ModelDiretory/" to "SuperDirectory/"
	where SuperDirectory is the directory that contains the master
	directory and all of your sim run directories.  Make sure the
	directory ends with a slash.

3) Use setup.q script to create multiple "sim" directories for your
simulation runs and load the appropriate model/stimulus
parameters. This will also copy your special file into the
directories.

4) Run the model using the command: ./special DBSstim.hoc quit.hoc
from inside each newly created directory.

5) DBSstim.hoc loads several parameter files at the beginning, which
is why you need the nk, vsh.txt, etc. files in each run directory.

6) If you would like to custumize the extracellular stimulus trains
played into the model, then you will need to generate Harmaline.dat,
Poisson.dat, and Stim_x.dat files in master directory using
generate_stim_poisson_harmaline.m file.
