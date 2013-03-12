From inet!cs.utexas.edu!cline Tue Oct 31 17:10:31 CST 1989
Received: from mojave.cs.utexas.edu by cs.utexas.edu (5.59/1.44)
	id AA29509; Tue, 31 Oct 89 17:11:51 CST
Posted-Date: Tue, 31 Oct 89 17:10:31 CST
Message-Id: <8910312310.AA04442@mojave.cs.utexas.edu>
Received: by mojave.cs.utexas.edu (14.5/1.4-Client)
	id AA04442; Tue, 31 Oct 89 17:10:34 cst
Date: Tue, 31 Oct 89 17:10:31 CST
X-Mailer: Mail User's Shell (6.5 4/17/89)
From: cline@cs.utexas.edu (Alan Cline)
To: ehg@research.att.com
Subject: New FITPACK Subset for netlib


This new version of FITPACK distributed by netlib is about 20% of 
the total package in terms of characters, lines of code, and num-
ber of subprograms. However, these 25 subprograms represent about
95% of usages of the package.  What has been omitted are such ca-
pabilities as:
  1. Automatic tension determination,
  2. Derivatives, arclengths, and enclosed areas for planar 
     curves,
  3. Three dimensional curves,
  4. Special surface fitting using equispacing assumptions,
  5. Surface fitting in annular, wedge, polar, toroidal, lunar,
     and spherical geometries,
  6. B-splines in tension generation and usage,
  7. General surface fitting in three dimensional space.

(The code previously circulated in netlib is less than 10% of the
total  package  and is more than a decade old.  Its usage is dis-
couraged.)

Please note:  Two versions of the subroutine snhcsh are included.
Both serve the same purpose:  obtaining approximations to certain
hyperbolic trigonometric-like functions.  The first is less accu-
rate (but more efficient) than the second.  Installers should se- 
lect the one with the precision they desire.

Interested parties can obtain the entire package on disk or  tape
from Pleasant  Valley Software, 8603 Altus Cove, Austin TX (USA),
78759 at a cost of $495 US. A 340 page manual  is  available  for
$30  US  per  copy.  The  package  includes  examples and machine
readable documentation.
