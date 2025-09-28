//
// C++ Implementation: About
//
//
// Author: Alexander Stein <webmaster@alexandersteinchanneler1978.com>, (C) 2025
//
// Copyright: See README file that comes with this source code
//
//

#include "about.h"
#include "constants.h"

#include <QPushButton>
#include <QGridLayout>
#include <QTextBrowser>

About::About()
{
    QGridLayout *L = new QGridLayout(this);
    browser = new QTextBrowser(this);
    browser->setOpenExternalLinks(true);
    L->addWidget(browser, 0, 0);
    L->addWidget(close = new QPushButton("Close", this), 1, 0);
    connect(close, SIGNAL(clicked()), this, SIGNAL(closeThis()));
    browser->setHtml(QString("\
<html>\
<head>\
<title>About MolSpektanalysis</title>\
</head>\
<body>\
<h1 align=\"center\">MolSpektanalysis ") + MAVersion + "</h1>\
<h6 align=\"center\">&copy; Alexander Stein</h6>\
<p>This is a program I have developed during the work on my diploma and PHD theses in the\
   group of Prof. Tiemann in the <a href=\"https://www.iqo.uni-hannover.de/\">Institute of Quantum Optics</a> as a tool for the analysis of my spectroscopic results.\
   The programs <i>ClassCluster</i> and <i>DrawSound</i> use significant parts of the internal libraries of <i>MolSpektAnalysis</i></p>\
<p>The program was partly developed under huge time pressure and is thus by far not free from bugs.\
   Because of this I do not provide any warranty for the use, if you use it you do it purely on your own risk.\
   I give you the advise to save your results regularly, this program is not free from crashes.\
   If intended to be used seriously, you should double check at least a significant sample of your results obtained with this program using\
   an independent method.</p>\
<p>The program is made using <a href=\"https://www.qt.io\">Qt</a>.\
   For numerical calculations, especially any kinds of fits, I often exercised <a href=\"http://numerical.recipes\">Numerical Recipes</a>.\
   To solve the molecular Schr&ouml;dingers equation (central field) the algorithm from\
   <a href=\"http://www.ams.org/journals/mcom/1961-15-076/S0025-5718-1961-0129566-X/\">J. W. Cooley, <i>Mathematics of Computation</i>, <b>15</b>, 363 (1961)</a>\
   is used.</p>\
<p>The program is free to use, it is licensed under the <a href=\"https://www.gnu.org/licenses/lgpl-3.0.de.html\">GNU LGPL v3</a> and the source code can be found on <a href=\"https://github.com/AlexanderStein1978/AlexanderStein1978\" Github</a>\
   If you have suggestions for improvements you can send it to me, but I cannot guarantee that I will find the time to implement all.\
   The same is true for bugs you observe. Especially, if it is a bug which significantly hampers your work, if it is a bug which leads to wrong results\
   or is a reproducible crash, I would be grateful for the information. The more significant information you provide, the easier you make it to me to\
   reproduce the bug, the more likely is it, that I will find the time to fix it and will provide you the improved version.</p>\
<p>I do not mind if you ask me questions how to use this program.</p>\
<p>If this program should significantly help you creating results you use for publications,\
   I would be grateful if you could cite one of our articles significant parts of the results published in were created under the use of this program:</p>\
   You can reach to me at <a href=\"mailto:webmaster@alexandersteinchanneler1978.com\"</a>.\
<ul>\
<li><a href=\"http://link.springer.com/article/10.1140/epjd%2Fe2008-00089-y\">A. Stein, A. Pashov, P. F. Staanum, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The European Physical Journal D</i>, <b>48</b>, 177 (2008)</a></li>\
<li><a href=\"http://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.042508\">A. Stein, H. Kn&ouml;ckel and E. Tiemann, <i>Physical Review A</i>, <b>78</b>\
    042508 (2008)</a></li>\
<li><a href=\"http://link.springer.com/article/10.1140/epjd%2Fe2010-00058-y\">A. Stein, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The European Physical Journal D</i>, <b>57</b>, 171 (2010)</a></li>\
<li><a href=\"http://link.springer.com/article/10.1140/epjd%2Fe2011-20229-6\">A. Stein, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The European Physical Journal D</i>, <b>64</b>, 227 (2011)</a></li>\
<li><a href=\"http://scitation.aip.org/content/aip/journal/jcp/134/2/10.1063/1.3524312\">M. Ivanova, A. Stein, A. Pashov, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The Journal of Chemical Physics</i>, <b>134</b>, 024321 (2011)</a></li>\
<li><a href=\"http://scitation.aip.org/content/aip/journal/jcp/135/17/10.1063/1.3652755\">\
    M. Ivanova, A. Stein, A. Pashov, A. V. Stolyarov, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The Journal of Chemical Physics</i>, <b>135</b>, 174303 (2011)</a></li>\
<li><a href=\"http://scitation.aip.org/content/aip/journal/jcp/138/9/10.1063/1.4793315\">M. Ivanova, A. Stein, A. Pashov, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The Journal of Chemical Physics</i>, <b>138</b>, 094315 (2013)</a></li>\
<li><a href=\"http://scitation.aip.org/content/aip/journal/jcp/138/11/10.1063/1.4795205\">A. Stein, M. Ivanova, A. Pashov, H. Kn&ouml;ckel and E. Tiemann,\
    <i>The Journal of Chemical Physics</i>, <b>138</b>, 114306 (2013)</a></li>\
</ul>\
</body>\
</html>");
}
