# SAMBA is a computer program that solves the problem of leakage between connected aquifers using a semi-analytic method
# Copyright (C) 2012 Arnaud REVEILLERE
# 
# Contact: a.reveillere at brgm.fr / arnaud.reveillere @ centraliens.net
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

from math import log

def float2(s):
    '''Convert a string s to a float. Deals with '0.3149-124' (=0.3149e-124) format.'''
    
    ss = ''
    for i in range(len(s)):
        if i>=1 and s[i] == '-' and s[i-1] not in ['e', 'E']:
            ss += 'e-'
        else:
            ss += s[i]
            
    return float(ss)

    
def GetBrineMassFlowFromCOFT(path, dt_min = 0):
    '''Return [[time(s)], [mass flow rate Connexion 1 (kg/s)], [mass flow rate Connexion 2]...]
    from the COFT file in path (relative or absolute) with a minimum time interval dt_min (s).
    This works with TOUGH2/eco2n output COFT files'''

    # Initialization
    l = 47 # COFT parameter
    time = 0
    lasttime = - dt_min
    # Read the Tough2 COFT output file
    f = open(path+'/COFT', 'r')
    lines = f.readlines()
    f.close()

    l = 47 # number of characters of a connexion in COFT
    NumberOfConnections = int((len(lines[0]) - 27)/73  + 1) # number of connections recorded in the COFT file
    conns = [[]]*(NumberOfConnections+1) # the +1 is for the time in the first column
    
    for line in lines:
        time = float2(line[8:20])
        
        if time - lasttime >= dt_min:
            lasttime = time
            conns[0] = conns[0]  + [time]
            for c in range(NumberOfConnections):
                conns[c+1] = conns[c+1] + [float2(line[25 + l*c + 19: 25 + l*c + 29])]
    return conns


    
def BuildExpRefinement(x, y, Shape, a_x, a_y, z_min, z_max, DivisionNumber, CentralArea, StripWidth):
    '''creates the input files of LGR.exe (for Local Grid Refinement) mesh-making scripts developed by Audigane et al. (2011)
    in order to build an logarithmic refinement in a TOUGH2 MESH file. Cf. input parameters definition in Reveillere, 2013.'''
    
    xn = Shape[0]
    xp = Shape[1]
    yn = Shape[2]
    yp = Shape[3]
        
    g = DivisionNumber
    
    n_xn = CentralArea[0]
    n_xp = CentralArea[1]
    n_yn = CentralArea[2]
    n_yp = CentralArea[3]

    l = StripWidth


    # Compute lateral Extensions
    Ext_xp = xp*l*a_x
    Ext_yp = yp*l*a_y
    Ext_xn = xn*l*a_x
    Ext_yn = yn*l*a_y

    # Compute offsets due to the internal area
    D_xp = n_xp*a_x
    D_xn = n_xn*a_x
    D_yp = n_yp*a_y
    D_yn = n_yn*a_y


    # Write the local grid refinement file

    f = open('LGR.itt','w')
    f.write("! nb of domain you want to remesh\n%s\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % g)
    f.write('%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, 2, 2, 1\n' % (x - D_xn - Ext_xn, x + D_xp + Ext_xp, y - D_yn - Ext_yn, y + D_yp + Ext_yp, z_min, z_max))

    for gg in range(1, g):
        f.write('%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, 2, 2, 1\n' % (x - D_xn - Ext_xn*2**-gg, x + D_xp + Ext_xp*2**-gg, y - D_yn - Ext_yn*2**-gg, y + D_yp + Ext_yp*2**-gg, z_min, z_max))

    f.close()

    #call("LGR.exe")

    
def BuildExponentialMesh(StripWidth, SymetryCase, LateralExtension, n, a_int, LayersDepths):
    '''creates the input files of Femesh2tough.exe and LGR.exe (for Local Grid Refinement) mesh-making 
    scripts developed by Audigane et al. (2011) in order to build an exponential TOUGH2 MESH file.
    Cf. input parameters definition in Reveillere, 2013.'''
    
    # Compute the internal area extent
    D = a_int*2**n
    
    # Compute the number g of LGRs that will be done
    g = int(log((LateralExtension-D)/(2*StripWidth*a_int))/log(2) + 1)
    
    # Size of the external cells (m)
    a_ext = a_int*2**g

    # minimal, maximal depth and number of layers
    zmin = str(LayersDepths[0])
    zmax = str(LayersDepths[-1])
    nk = len(LayersDepths)-1
    
    # Prepare the text for the depth discretization
    depth_discr = "\npoint k\n"
    for z in reversed(LayersDepths):
        depth_discr += str(z)+'\n'
        
    # Compute the real extension (useless, it's just for information... :-))
    Ext_real = a_int * 2*StripWidth*2**g + D
    
    # Case 1: x>=0; y>=0

    if SymetryCase == 1:
            
        if n != 0:

            # First mesh
            f = open('Femesh2tough.itt','w')

            if n <= g : # D <= a_ext

                f.write('ni\n')
                f.write('%d\n' % (2*StripWidth + 1,))

                f.write('nj\n')
                f.write('%d\n' % (2*StripWidth + 1,))

                f.write('nk\n')
                f.write('%d\n' % (nk,))
                
                    
                f.write("\npoint i\n")
                f.write('%.1f\n' % 0)
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext+ D,))
              
                f.write("\npoint j\n")
                f.write('%.1f\n' % 0)
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext+ D,))

                f.write(depth_discr)

            elif n > g : # D > a_ext

                f.write('ni\n')
                f.write('%d\n' % (2*StripWidth + 1 + 2**(n-g),))

                f.write('nj\n')
                f.write('%d\n' % (2*StripWidth + 1 + 2**(n-g),))

                f.write('nk\n')
                f.write('%d\n' % (nk,))
                
                f.write("\npoint i\n")
                for i in range(0, 2**(n-g) + 1):
                    f.write('%.1f\n' % (i*a_ext,))              
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext+ D,))
              
                f.write("\npoint j\n")
                for i in range(0, 2**(n-g) + 1):
                    f.write('%.1f\n' % (i*a_ext,))   
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write(depth_discr)

            f.close()
            #call("Femesh2tough.exe")            
                
            

            # LGRs

            gig_file =""
            nb_gig = 0

            for gg in range(1, g+1):

                r = a_int*StripWidth*2**(g+1-gg)

                if gg <= g-n :
                    gig_file += "%.1f, %s, %s, %.1f, %s, %s, 2, 1, 1\n" % (D, D + r, D, zmin, zmax) # zone 2
                    gig_file += "0., %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (D, D, D + r, zmin, zmax) # zone 3
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, D, D + r, zmin, zmax) # zone 5
                    nb_gig += 3
                    
                else:
                    gig_file += '0., %s, %s, %.1f, %s, %s, 2, 2, 1\n' % (D + r, D + r, zmin, zmax)
                    nb_gig += 1

            gig_file = "! nb of domain you want to remesh\n%d\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % nb_gig + gig_file


            f = open('LGR.itt','w')
            f.write(gig_file)
            f.close()

            #call("LGR.exe")

            

        elif n == 0:


            # Number of cells side to side
            Nm = nk * 2*StripWidth*2**g

            
            # First mesh

            f = open('Femesh2tough.itt','w')

            f.write('ni\n')
            f.write('%d\n' % (2*StripWidth,))

            f.write('nj\n')
            f.write('%d\n' % (2*StripWidth,))

            f.write('nk\n')
            f.write('%d\n' % (nk,))


            f.write("\npoint i\n")
            for i in range(0, 2*StripWidth+1):
                f.write('%.1f\n' % (i*a_ext,))

            f.write("\npoint j\n")
            for i in range(0, 2*StripWidth+1):
                f.write('%.1f\n' % (i*a_ext,))

            f.write(depth_discr)

            f.close()

            #call("Femesh2tough.exe")
            
            # LGRs
            f = open('LGR.itt','w')

            f.write("! nb of domain you want to remesh\n%d\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % g)

            for gg in range(1, g+1):
                r = a_int*StripWidth*2**(g+1-gg)
                f.write("0., %s, %s, %.1f, %s, %s, 2, 2, 1\n" % (r, r, zmin, zmax))

            f.close()

            #call("LGR.exe")





    # Case 2: y>=0

    elif SymetryCase == 2:

        if n != 0:

            # First Mesh
            f = open('Femesh2tough.itt','w')
            
            if n < g : # D < a_ext

                f.write('ni\n')
                f.write('%d\n' % (4*StripWidth + 1,))

                f.write('nj\n')
                f.write('%d\n' % (2*StripWidth + 1,))

                f.write('nk\n')
                f.write('%d\n' % (nk,))


                f.write("\npoint i\n")
                for i in range(-2*StripWidth, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write("\npoint j\n")
                f.write('%.1f\n' % 0)
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write(depth_discr)

           
            elif n == g : # D = a_ext

                f.write('ni\n')
                f.write('%d\n' % (4*StripWidth + 2,))

                f.write('nj\n')
                f.write('%d\n' % (2*StripWidth + 1,))

                f.write('nk\n')
                f.write('%d\n' % (nk,))


                f.write("\npoint i\n")
                for i in range(-2*StripWidth, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                f.write('.0\n')
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write("\npoint j\n")
                f.write('%.1f\n' % 0)
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext+ D,))

                f.write(depth_discr)

      
            elif n > g : # D > a_ext

                f.write('ni\n')
                f.write('%d\n' % (4*StripWidth + 2 + 2**(n+1-g),))

                f.write('nj\n')
                f.write('%d\n' % (2*StripWidth + 1 + 2**(n-g),))

                f.write('nk\n')
                f.write('%d\n' % (nk,))


                f.write("\npoint i\n")
                for i in range(-2*StripWidth, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                for i in range(-2**(n-g), 2**(n-g)+1):
                    f.write('%.1f\n' % (i*a_ext,))            
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write("\npoint j\n")
                for i in range(0, 2**(n-g) + 1):
                    f.write('%.1f\n' % (i*a_ext,))  
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext+ D,))

                f.write(depth_discr)                 
            f.close()
            #call("Femesh2tough.exe")

                   

            # LGRs

            gig_file =""
            nb_gig = 0

            for gg in range(1, g+1):

                r = a_int*StripWidth*2**(g+1-gg)

                if gg < g-n :
                    gig_file += "%.1f, %.1f, 0., %.1f, %s, %s, 2, 1, 1\n" % (-r - D, -D, D, zmin, zmax) # zone 1
                    gig_file += "%.1f, %.1f, 0., %.1f, %s, %s, 2, 1, 1\n" % (D, r + D, D, zmin, zmax) # zone 2
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (-D, D, D, r + D, zmin, zmax) # zone 3
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r - D, -D, D, r + D, zmin, zmax) # zone 4
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, D, r + D, zmin, zmax) # zone 5
                    nb_gig += 5
                    
                elif gg == g-n:
                    gig_file += "%.1f, %.1f, 0., %.1f, %s, %s, 2, 1, 1\n" % (-r - D, -D, D, zmin, zmax) # zone 1
                    gig_file += "%.1f, %.1f, 0., %.1f, %s, %s, 2, 1, 1\n" % (D, r + D, D, zmin, zmax) # zone 2
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (-D, D, D, r + D, zmin, zmax) # zone 3
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r - D, -D, D, r + D, zmin, zmax) # zone 4
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, D, r + D, zmin, zmax) # zone 5
                    gig_file += "%.1f, %.1f, 0., %.1f, %s, %s, 2, 1, 1\n" % (-D, D, r + D, zmin, zmax) # zone 3+6
                    nb_gig += 6
                    
                else:
                    gig_file += '%.1f, %.1f, 0., %.1f, %s, %s, 2, 2, 1\n' % (-r - D, r + D, r + D, zmin, zmax)
                    nb_gig += 1

            gig_file = "! nb of domain you want to remesh\n%d\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % nb_gig + gig_file


            f = open('LGR.itt','w')
            f.write(gig_file)
            f.close()

            #call("LGR.exe")

            

        elif n == 0:
            
            # Affiche le nb de mailles final hors extension
            Nm = nk * 4*StripWidth*2**g
            print('%d cells' % Nm)
            
            #  First mesh

            f = open('Femesh2tough.itt','w')

            f.write('ni\n')
            f.write('%d\n' % (4*StripWidth,))

            f.write('nj\n')
            f.write('%d\n' % (2*StripWidth,))

            f.write('nk\n')
            f.write('%d\n' % (nk,))


            f.write("\npoint i\n")
            for i in range(-2*StripWidth, 2*StripWidth+1):
                f.write('%.1f\n' % (i*a_ext,))

            f.write("\npoint j\n")
            for i in range(0, 2*StripWidth+1):
                f.write('%.1f\n' % (i*a_ext,))

            f.write(depth_discr)

            f.close()

            #call("Femesh2tough.exe")
            
            # LGRs
            f = open('LGR.itt','w')

            f.write("! nb of domain you want to remesh\n%d\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % g)

            for gg in range(1, g+1):
                
                r = a_int*StripWidth*2**(g+1-gg)
                f.write("%.1f, %.1f, 0., %.1f, %s, %s, 2, 2, 1\n" % (-r, r, r, zmin, zmax))

            f.close()

            #call("LGR.exe")



    # Case 3: full mesh
    elif SymetryCase == 3:

        if n != 0:

            # First Mesh
            f = open('Femesh2tough.itt','w')
            
            if n < g : # D < a_ext

                f.write('ni\n')
                f.write('%d\n' % (4*StripWidth + 1,))

                f.write('nj\n')
                f.write('%d\n' % (4*StripWidth + 1,))

                f.write('nk\n')
                f.write('%d\n' % (nk,))


                f.write("\npoint i\n")
                for i in range(-2*StripWidth, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write("\npoint j\n")
                for i in range(-2*StripWidth, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                for i in range(0, 2*StripWidth+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write(depth_discr)

      
            elif n >= g : # D >= a_ext

                f.write('ni\n')
                f.write('%d\n' % (4*l + 2 + 2**(n+1-g),))

                f.write('nj\n')
                f.write('%d\n' % (4*l + 2 + 2**(n+1-g),))

                f.write('nk\n')
                f.write('%d\n' % (nk,))


                f.write("\npoint i\n")
                for i in range(-2*l, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                for i in range(-2**(n-g), 2**(n-g)+1):
                    f.write('%.1f\n' % (i*a_ext,))            
                for i in range(0, 2*l+1):
                    f.write('%.1f\n' % (i*a_ext + D,))
                    
                f.write("\npoint j\n")
                for i in range(-2*l, 1):
                    f.write('%.1f\n' % (i*a_ext - D,))
                for i in range(-2**(n-g), 2**(n-g)+1):
                    f.write('%.1f\n' % (i*a_ext,))            
                for i in range(0, 2*l+1):
                    f.write('%.1f\n' % (i*a_ext + D,))

                f.write(depth_discr)
                    
            f.close()
            #call("Femesh2tough.exe")

                   

            # LGRs

            gig_file =""
            nb_gig = 0

            for gg in range(1, g+1):

                r = a_int*l*2**(g+1-gg)

                if gg < g-n :
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 1, 1\n" % (-r - D, -D, -D, D, zmin, zmax) # zone 1
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 1, 1\n" % (D, r + D, -D, D, zmin, zmax) # zone 2
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (-D, D, D, r + D, zmin, zmax) # zone 3
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r - D, -D, D, r + D, zmin, zmax) # zone 4
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, D, r + D, zmin, zmax) # zone 5
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r - D, -D, -D, -r - D, zmin, zmax) # zone 7
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (-D, D, -D, -r - D, zmin, zmax) # zone 8
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, -D, -r - D, zmin, zmax) # zone 9
                    nb_gig += 8
                    
                elif gg == g-n:
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 1, 1\n" % (-r - D, -D, -D, D, zmin, zmax) # zone 1
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 1, 1\n" % (D, r + D, -D, D, zmin, zmax) # zone 2
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (-D, D, D, r + D, zmin, zmax) # zone 3
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r - D, -D, D, r + D, zmin, zmax) # zone 4
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, D, r + D, zmin, zmax) # zone 5
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-D, D, -D, D, zmin, zmax) # zone 6
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r - D, -D, -D, -r - D, zmin, zmax) # zone 7
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 1, 2, 1\n" % (-D, D, -D, -r - D, zmin, zmax) # zone 8
                    gig_file += "%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (D, r + D, -D, -r - D, zmin, zmax) # zone 9
                    nb_gig += 9
                    
                else:
                    gig_file += '%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n' % (-r - D, r + D, -r - D, r + D, zmin, zmax)
                    nb_gig += 1

            gig_file = "! nb of domain you want to remesh\n%d\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % nb_gig + gig_file


            f = open('LGR.itt','w')
            f.write(gig_file)
            f.close()

            #call("LGR.exe")

            

        elif n == 0:
            
            # Affiche le nb de mailles final hors extension
            Nm = nk * 8*l*2**g
            print('%d cells' % Nm)
            
            # First mesh

            f = open('Femesh2tough.itt','w')

            f.write('ni\n')
            f.write('%d\n' % (4*l,))

            f.write('nj\n')
            f.write('%d\n' % (4*l,))

            f.write('nk\n')
            f.write('%d\n' % (nk,))


            f.write("\npoint i\n")
            for i in range(-2*l, 2*l+1):
                f.write('%.1f\n' % (i*a_ext,))

            f.write("\npoint j\n")
            for i in range(-2*l, 2*l+1):
                f.write('%.1f\n' % (i*a_ext,))

            f.write(depth_discr)

            f.close()

            #call("Femesh2tough.exe")
            
            # LGRs
            f = open('LGR.itt','w')

            f.write("! nb of domain you want to remesh\n%d\n! xmin, xmax, ymin, ymax, zmin, zmax, xfact (integer>=1), yfact (integer>=1), zfact (integer>=1)\n" % g)

            for gg in range(1, g+1):
                
                r = a_int*l*2**(g+1-gg)
                f.write("%.1f, %.1f, %.1f, %.1f, %s, %s, 2, 2, 1\n" % (-r, r, -r, r, zmin, zmax))

            f.close()

            #call("LGR.exe")

			
"""
References:

Audigane, P., Chiaberge, C., Mathurin, F., Lions, J., Picot-Colbeaux, G.: A work?ow for handling heterogeneous
3D models with the TOUGH2family of codes: applications to numerical modeling ofCO2 geological
storage. Comput. Geosci. 37, 610–662 (2011)

Reveillere, A., 2013. SAMBA v1.0 – User Guide. BRGM public report, in preparation.
"""
