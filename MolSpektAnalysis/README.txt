About this program:

MolSpektAnalysis 

Author: Alexander Stein

This is a program I have developed during the work on my diploma and PHD thesis in the
group of Prof. Tiemann as a tool for the analysis of my spectroscopic results.

The program was partly developed under huge time pressure and is thus by far not free from bugs.
Because of this I do not provide any warranty for the use, if you use it you do it purely on your own risk.

I give you the advise to save your results regularly, this program is not free from crashes.

If intended to be used seriously, you should double check at least a significant sample of your results
obtained with this program using an independent method.

The program is made using Qt.
For numerical calculations, especially any kinds of fits, I often exercised Numerical Recipes.
To solve the molecular Schrödingers equation (central field) the algorithm from
J. W. Cooley, Mathematics of Computation, 15, 363 (1961)
is used.
The program is free to use, it is licensed under the GNU LGPL v3 which can be found at "https://www.gnu.org/licenses/lgpl-3.0.de.html">. The source code is to be found at "https://github.com/AlexanderStein1978/AlexanderStein1978" 
If you have suggestions for improvements you can send it to me, but I cannot guarantee that I will find the time to implement all.
The same is true for bugs you observe. Especially, if it is a bug which significantly hampers your work, if it is a bug which leads to wrong results
or is a reproducible crash, I would be grateful for the information. The more significant information you provide, the easier you make it to me to
reproduce the bug, the more likely is it, that I will find the time to fix it and will provide you the improved version.
I do not mind if you ask me questions how to use this program.
If this program should significantly help you creating results you use for publications,
I would be grateful if you could cite one of our articles significant parts of the results published in were created under the use of this program:
You can reach to me at "webmaster@alexandersteinchanneler1978.com"

I do not mind if you ask me questions how to use this program.

If this program should significantly help you creating results you use for publications,
I would be grateful if you could cite one of our articles significant parts of the results published in were
created under the use of this program:

A. Stein, A. Pashov, P. F. Staanum, H. Knöckel and E. Tiemann, The European Physical Journal D, 48, 177 (2008)
A. Stein, H. Knöckel and E. Tiemann, Physical Review A, 78, 042508 (2008)
A. Stein, H. Knöckel and E. Tiemann, The European Physical Journal D, 57, 171 (2010)
A. Stein, H. Knöckel and E. Tiemann, The European Physical Journal D, 64, 227 (2011)
M. Ivanova, A. Stein, A. Pashov, H. Knöckel and E. Tiemann, The Journal of Chemical Physics, 134, 024321 (2011)
M. Ivanova, A. Stein, A. Pashov, A. V. Stolyarov, H. Knöckel and E. Tiemann, The Journal of Chemical Physics, 135, 174303 (2011)
M. Ivanova, A. Stein, A. Pashov, H. Knöckel and E. Tiemann, The Journal of Chemical Physics, 138, 094315 (2013)
A. Stein, M. Ivanova, A. Pashov, H. Knöckel and E. Tiemann, The Journal of Chemical Physics, 138, 114306 (2013)


To the source code:

Please excuse the partly quite bad readability. Because of the already mentioned time pressure and the fact, that I had at
the beginning of the developement of this program a lack of experience in reading code from other authors or even reading
the source code of this program after years again.

You are allowed to modify the source code, but please keep the current text of this file and the html code in the file
about.cpp and add to each version of this text a note about your modifications if you make your modified version available
to third persons.

Additionally, I would be grateful if you could send me a copy of the files you modified. This is a way you could
reciprocate a little that I made this program available to you. If you do not write that you do not want this and
if I find that your code is an improvement I would add it to my version of this program and also mention you with the
name you tell me for this purpose in the 'about' text. This way you could also benefit from my improvements or the
improvements third persons send me without having to merge your own improvements again and again.

To be able to compile the current version of this code, you need the freely available Qt library (tested are versions
between 4.8 and 5.6 on Linux and Windows platforms) and the commercially available source code package of Numerical Recipes Third Edition. 
If you do not want to purchase the commercial library, you have the choice to comment out the code using it (in fit.h
and fit.cpp) removing the possibility to fit Dunham coefficients and potentials, replace it by alternatives or, again,
send your changes to me. If the code changed by you directly compiles without any problems or if I should find the time
fix the bugs hampering it I will send you back the Linux or Windows binary you asked me for.   
  