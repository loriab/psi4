"""
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.


- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'<subset>'`` <members_description>

"""
import re
import qcdb

# <<< DARC Database Module >>>
dbse = 'DARC'

# <<< Database Members >>>
HRXN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', ]
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, '1'                     )] = ['%s-%s-reagent'      % (dbse, 'ethene'),
                                                               '%s-%s-reagent'      % (dbse, 'butadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'P1') ]
RXNM['%s-%s'            % (dbse, '1'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '1')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '2'                     )] = ['%s-%s-reagent'      % (dbse, 'ethine'),
                                                               '%s-%s-reagent'      % (dbse, 'butadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'P2') ]
RXNM['%s-%s'            % (dbse, '2'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '2')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '3'                     )] = ['%s-%s-reagent'      % (dbse, 'ethene'),
                                                               '%s-%s-reagent'      % (dbse, 'cyclopentadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'P3') ]
RXNM['%s-%s'            % (dbse, '3'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '3')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '4'                     )] = ['%s-%s-reagent'      % (dbse, 'ethine'),
                                                               '%s-%s-reagent'      % (dbse, 'cyclopentadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'P4') ]
RXNM['%s-%s'            % (dbse, '4'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '4')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '5'                     )] = ['%s-%s-reagent'      % (dbse, 'ethene'),
                                                               '%s-%s-reagent'      % (dbse, 'cyclohexadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'P5') ]
RXNM['%s-%s'            % (dbse, '5'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '5')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '6'                     )] = ['%s-%s-reagent'      % (dbse, 'ethine'),
                                                               '%s-%s-reagent'      % (dbse, 'cyclohexadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'P6') ]
RXNM['%s-%s'            % (dbse, '6'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '6')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '7'                     )] = ['%s-%s-reagent'      % (dbse, 'furane'),
                                                               '%s-%s-reagent'      % (dbse, 'maleine'),
                                                               '%s-%s-reagent'      % (dbse, 'P7endo') ]
RXNM['%s-%s'            % (dbse, '7'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '7')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '8'                     )] = ['%s-%s-reagent'      % (dbse, 'furane'),
                                                               '%s-%s-reagent'      % (dbse, 'maleine'),
                                                               '%s-%s-reagent'      % (dbse, 'P7exo') ]
RXNM['%s-%s'            % (dbse, '8'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '8')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '9'                     )] = ['%s-%s-reagent'      % (dbse, 'furane'),
                                                               '%s-%s-reagent'      % (dbse, 'maleinimide'),
                                                               '%s-%s-reagent'      % (dbse, 'P8endo') ]
RXNM['%s-%s'            % (dbse, '9'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '9')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '10'                    )] = ['%s-%s-reagent'      % (dbse, 'furane'),
                                                               '%s-%s-reagent'      % (dbse, 'maleinimide'),
                                                               '%s-%s-reagent'      % (dbse, 'P8exo') ]
RXNM['%s-%s'            % (dbse, '10'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '10')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '11'                    )] = ['%s-%s-reagent'      % (dbse, 'cyclopentadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'maleine'),
                                                               '%s-%s-reagent'      % (dbse, 'P9endo') ]
RXNM['%s-%s'            % (dbse, '11'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '11')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '12'                    )] = ['%s-%s-reagent'      % (dbse, 'cyclopentadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'maleine'),
                                                               '%s-%s-reagent'      % (dbse, 'P9exo') ]
RXNM['%s-%s'            % (dbse, '12'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '12')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '13'                    )] = ['%s-%s-reagent'      % (dbse, 'cyclopentadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'maleinimide'),
                                                               '%s-%s-reagent'      % (dbse, 'P10endo') ]
RXNM['%s-%s'            % (dbse, '13'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '13')], [-1, -1, 1]))

ACTV['%s-%s'            % (dbse, '14'                    )] = ['%s-%s-reagent'      % (dbse, 'cyclopentadiene'),
                                                               '%s-%s-reagent'      % (dbse, 'maleinimide'),
                                                               '%s-%s-reagent'      % (dbse, 'P10exo') ]
RXNM['%s-%s'            % (dbse, '14'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '14')], [-1, -1, 1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1'                     )] =   -43.8
BIND['%s-%s'            % (dbse, '2'                     )] =   -59.3 
BIND['%s-%s'            % (dbse, '3'                     )] =   -30.0 
BIND['%s-%s'            % (dbse, '4'                     )] =   -33.1 
BIND['%s-%s'            % (dbse, '5'                     )] =   -36.5 
BIND['%s-%s'            % (dbse, '6'                     )] =   -48.2 
BIND['%s-%s'            % (dbse, '7'                     )] =   -14.4
BIND['%s-%s'            % (dbse, '8'                     )] =   -16.2 
BIND['%s-%s'            % (dbse, '9'                     )] =   -17.2
BIND['%s-%s'            % (dbse, '10'                    )] =   -19.2
BIND['%s-%s'            % (dbse, '11'                    )] =   -31.6
BIND['%s-%s'            % (dbse, '12'                    )] =   -32.1
BIND['%s-%s'            % (dbse, '13'                    )] =   -34.1
BIND['%s-%s'            % (dbse, '14'                    )] =   -34.4
#Taken from Johnson, E. R.; Mori-Sanchez, P.; Cohen, A. J; Yang, W. J. Chem. Phys. 2008, 129, 204112.
#(estimated CCSD(T)/CBS); all values are in kcal/mol.
# collected into dbse in 10.1021/ct900489g Table S17

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Diels-Alder: ethene + butadiene -->"""
TAGL['%s-%s'            % (dbse, '2'                     )] = """Diels-Alder: ethine + butadiene -->"""
TAGL['%s-%s'            % (dbse, '3'                     )] = """Diels-Alder: ethene + cyclopentadiene -->"""
TAGL['%s-%s'            % (dbse, '4'                     )] = """Diels-Alder: ethine + cyclopentadiene -->"""
TAGL['%s-%s'            % (dbse, '5'                     )] = """Diels-Alder: ethene + cyclohexadiene -->"""
TAGL['%s-%s'            % (dbse, '6'                     )] = """Diels-Alder: ethine + cyclohexadiene -->"""
TAGL['%s-%s'            % (dbse, '7'                     )] = """Diels-Alder: furane + maleine --> (endo)"""
TAGL['%s-%s'            % (dbse, '8'                     )] = """Diels-Alder: furane + maleine --> (exo)"""
TAGL['%s-%s'            % (dbse, '9'                     )] = """Diels-Alder: furane + maleinimide --> (endo)"""
TAGL['%s-%s'            % (dbse, '10'                    )] = """Diels-Alder: furane + maleinimide --> (exo)"""
TAGL['%s-%s'            % (dbse, '11'                    )] = """Diels-Alder: cyclopentadiene + maleine --> (endo)"""
TAGL['%s-%s'            % (dbse, '12'                    )] = """Diels-Alder: cyclopentadiene + maleine --> (exo)"""
TAGL['%s-%s'            % (dbse, '13'                    )] = """Diels-Alder: cyclopentadiene + maleinimide --> (endo)"""
TAGL['%s-%s'            % (dbse, '14'                    )] = """Diels-Alder: cyclopentadiene + maleinimide --> (exo)"""
TAGL['%s-%s-reagent'    % (dbse, 'P1'                    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P10endo'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P10exo'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P2'                    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P3'                    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P4'                    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P5'                    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P6'                    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P7endo'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P7exo'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P8endo'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P8exo'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P9endo'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'P9exo'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'butadiene'             )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cyclohexadiene'        )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cyclopentadiene'       )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ethene'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ethine'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'furane'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'maleine'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'maleinimide'           )] = """ """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'P1', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -1.250799702386    -0.112835135019    -2.634108874863
    C               -2.808528394402    -0.214575357044    -0.273479746780
    C               -1.306986864738     0.600056944295     2.044165853125
    C                1.306986864738    -0.600056944295     2.044165853125
    C                2.808528394402     0.214575357044    -0.273479746780
    C                1.250799702386     0.112835135019    -2.634108874863
    H               -2.249295767194    -0.224495244259    -4.425287485096
    H               -4.484801132739     0.979076443748    -0.498097811139
    H               -3.536213607854    -2.141567592374    -0.013824334102
    H               -2.334832408016     0.115037020126     3.769945251085
    H               -1.103163590900     2.661380918425     2.030687147770
    H                2.334832408016    -0.115037020126     3.769945251085
    H                1.103163590900    -2.661380918425     2.030687147770
    H                4.484801132739    -0.979076443748    -0.498097811139
    H                3.536213607854     2.141567592374    -0.013824334102
    H                2.249295767194     0.224495244259    -4.425287485096

""")

GEOS['%s-%s-%s' % (dbse, 'P10endo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.840681766982     0.165241552635    -2.123808730331
    C               -0.695382569347     1.519923502267    -1.452463672843
    C               -0.695382569347     1.519923502267     1.452463672843
    C                1.840681766982     0.165241552635     2.123808730331
    C                1.588312206006    -2.544184085500    -1.259342365738
    C                1.588312206006    -2.544184085500     1.259342365738
    C                3.509254791583     1.230474401835     0.000000000000
    C               -3.064166247652     0.110091940359    -2.207400300383
    C               -3.064166247652     0.110091940359     2.207400300383
    N               -4.285148076038    -0.571902336207     0.000000000000
    O               -3.802643887996    -0.368022718078     4.302030316661
    O               -3.802643887996    -0.368022718078    -4.302030316661
    H               -5.897662821822    -1.591488066078     0.000000000000
    H                2.446068934906     0.458466550545    -4.065968108972
    H               -0.745097889389     3.405211108406    -2.287972866738
    H               -0.745097889389     3.405211108406     2.287972866738
    H                2.446068934906     0.458466550545     4.065968108972
    H                1.255394083251    -4.136147910942    -2.498980065841
    H                1.255394083251    -4.136147910942     2.498980065841
    H                3.623335996454     3.294121126538     0.000000000000
    H                5.403887316301     0.417634994528     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'P10exo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.520158125453     0.553737673396    -2.121098779241
    C               -0.668213289700    -1.305177318946    -1.454201129240
    C               -0.668213289700    -1.305177318946     1.454201129240
    C                1.520158125453     0.553737673396     2.121098779241
    C                3.900486157755    -0.774166872023    -1.258793686515
    C                3.900486157755    -0.774166872023     1.258793686515
    C                1.232911700396     2.510023371332     0.000000000000
    C               -3.249676035478    -0.317688674781    -2.203332534848
    C               -3.249676035478    -0.317688674781     2.203332534848
    N               -4.574380039657     0.154967904341     0.000000000000
    O               -4.051780199159     0.027901317170     4.299756347188
    O               -4.051780199159     0.027901317170    -4.299756347188
    H               -6.344707410049     0.867112445017     0.000000000000
    H                1.470573037335     1.222936322564    -4.064627819680
    H               -0.390908414362    -3.149220753515    -2.330357442776
    H               -0.390908414362    -3.149220753515     2.330357442776
    H                1.470573037335     1.222936322564     4.064627819680
    H                5.224847632128    -1.716560898922    -2.500411943628
    H                5.224847632128    -1.716560898922     2.500411943628
    H               -0.584579250252     3.491931136907     0.000000000000
    H                2.759780971617     3.892443552518     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'P2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.352707270949     0.000000000000    -1.253692993531
    C                0.000000000000     0.000000000000    -2.814748982953
    C               -2.352707270949    -0.000000000000    -1.253692993531
    C               -2.352707270949    -0.000000000000     1.253692993531
    C                0.000000000000     0.000000000000     2.814748982953
    C                2.352707270949     0.000000000000     1.253692993531
    H                4.140341956918     0.000000000000    -2.263053165051
    H               -0.000000000000     1.640550343022    -4.089655961606
    H                0.000000000000    -1.640550343022    -4.089655961606
    H               -4.140341956918    -0.000000000000    -2.263053165051
    H               -4.140341956918    -0.000000000000     2.263053165051
    H                0.000000000000    -1.640550343022     4.089655961606
    H               -0.000000000000     1.640550343022     4.089655961606
    H                4.140341956918     0.000000000000     2.263053165051

""")

GEOS['%s-%s-%s' % (dbse, 'P3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -0.600184775679     0.382495924992    -2.116313437692
    C                0.974864496526     2.603721685225    -1.259286847739
    C                0.974864496526     2.603721685225     1.259286847739
    C               -0.600184775679     0.382495924992     2.116313437692
    C               -2.579565413726     0.302255567695     0.000000000000
    C                0.970472132130    -2.012439073393    -1.462135994565
    C                0.970472132130    -2.012439073393     1.462135994565
    H               -1.291015209751     0.445263624288    -4.054744880463
    H                2.083250126161     3.800771702678    -2.495039083125
    H                2.083250126161     3.800771702678     2.495039083125
    H               -1.291015209751     0.445263624288     4.054744880463
    H               -3.818434547538     1.951064436929     0.000000000000
    H               -3.705368612022    -1.431611494142     0.000000000000
    H                0.047696403531    -3.702260758932    -2.209182858285
    H                2.866601113725    -1.928407360099    -2.266593336914
    H                0.047696403531    -3.702260758932     2.209182858285
    H                2.866601113725    -1.928407360099     2.266593336914

""")

GEOS['%s-%s-%s' % (dbse, 'P4', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                0.000000000000    -2.106898586907    -0.401646694391
    C                2.330323873439    -1.254796423359     1.089419099179
    C                2.330323873439     1.254796423359     1.089419099179
    C               -0.000000000000     2.106898586907    -0.401646694391
    C                0.000000000000     0.000000000000    -2.428547639388
    C               -2.330323873439    -1.254796423359     1.089419099179
    C               -2.330323873439     1.254796423359     1.089419099179
    H                0.000000000000    -4.060288678517    -1.046040170379
    H                3.630776462205    -2.512871245543     2.039557542472
    H                3.630776462205     2.512871245543     2.039557542472
    H               -0.000000000000     4.060288678517    -1.046040170379
    H               -1.701327798261    -0.000000000000    -3.595992598839
    H                1.701327798261     0.000000000000    -3.595992598839
    H               -3.630776462205    -2.512871245543     2.039557542472
    H               -3.630776462205     2.512871245543     2.039557542472

""")

GEOS['%s-%s-%s' % (dbse, 'P5', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -0.000000000000     2.420205188102    -0.364967390978
    C               -2.349726910169     1.456175192680     1.058191551635
    C               -2.349726910169    -1.456175192680     1.058191551635
    C                0.000000000000    -2.420205188102    -0.364967390978
    C               -0.000000000000     1.257575647621    -2.949770510175
    C                0.000000000000    -1.257575647621    -2.949770510175
    C                2.349726910169     1.456175192680     1.058191551635
    C                2.349726910169    -1.456175192680     1.058191551635
    H               -0.000000000000     4.482067102726    -0.440332680189
    H               -2.322084184226     2.195634291903     2.989134565009
    H               -4.054787849674     2.190754610315     0.156971928312
    H               -4.054787849674    -2.190754610315     0.156971928312
    H               -2.322084184226    -2.195634291903     2.989134565009
    H                0.000000000000    -4.482067102726    -0.440332680189
    H               -0.000000000000     2.396878783958    -4.653525508571
    H                0.000000000000    -2.396878783958    -4.653525508571
    H                2.322084184226     2.195634291903     2.989134565009
    H                4.054787849674     2.190754610315     0.156971928312
    H                2.322084184226    -2.195634291903     2.989134565009
    H                4.054787849674    -2.190754610315     0.156971928312

""")

GEOS['%s-%s-%s' % (dbse, 'P6', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -0.000000000000     2.405178634118     0.455694309154
    C               -2.320477326687     1.254148714888     1.643291553505
    C               -2.320477326687    -1.254148714888     1.643291553505
    C                0.000000000000    -2.405178634118     0.455694309154
    C               -0.000000000000     1.451237392835    -2.340633098883
    C                0.000000000000    -1.451237392835    -2.340633098883
    C                2.320477326687     1.254148714888     1.643291553505
    C                2.320477326687    -1.254148714888     1.643291553505
    H               -0.000000000000     4.462841192983     0.531156614615
    H               -3.855059031457     2.416269170446     2.340890866279
    H               -3.855059031457    -2.416269170446     2.340890866279
    H                0.000000000000    -4.462841192983     0.531156614615
    H                1.664670746489     2.194817304974    -3.307291332227
    H               -1.664670746489     2.194817304974    -3.307291332227
    H                1.664670746489    -2.194817304974    -3.307291332227
    H               -1.664670746489    -2.194817304974    -3.307291332227
    H                3.855059031457     2.416269170446     2.340890866279
    H                3.855059031457    -2.416269170446     2.340890866279

""")

GEOS['%s-%s-%s' % (dbse, 'P7endo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.031333379610     0.277767902581    -2.007049485344
    C               -0.538070712552     1.656148232937    -1.438424924820
    C               -0.538070712552     1.656148232937     1.438424924820
    C                2.031333379610     0.277767902581     2.007049485344
    C                1.708818196932    -2.459579354561    -1.255790949302
    C                1.708818196932    -2.459579354561     1.255790949302
    O                3.518579498800     1.256918762520     0.000000000000
    C               -2.885711699719     0.223483903382    -2.156861174203
    C               -2.885711699719     0.223483903382     2.156861174203
    O               -4.129926564098    -0.544035046102     0.000000000000
    O               -3.675026940436    -0.264435851136     4.198424392842
    O               -3.675026940436    -0.264435851136    -4.198424392842
    H                2.858596872599     0.698353585306    -3.837654407512
    H               -0.551181179115     3.522513876290    -2.306430267009
    H               -0.551181179115     3.522513876290     2.306430267009
    H                2.858596872599     0.698353585306     3.837654407512
    H                1.356915615330    -4.010694153009    -2.536467443671
    H                1.356915615330    -4.010694153009     2.536467443671

""")

GEOS['%s-%s-%s' % (dbse, 'P7exo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.267390442973     0.995516873178    -2.006788937567
    C               -0.886054200589    -0.948064283336    -1.443067458502
    C               -0.886054200589    -0.948064283336     1.443067458502
    C                1.267390442973     0.995516873178     2.006788937567
    C                3.703153607785    -0.298023630974    -1.255017560477
    C                3.703153607785    -0.298023630974     1.255017560477
    O                0.894672417287     2.738133139260     0.000000000000
    C               -3.438982972294     0.107142897580    -2.147734362345
    C               -3.438982972294     0.107142897580     2.147734362345
    O               -4.782852461858     0.685599671759     0.000000000000
    O               -4.286061502245     0.450668011111     4.194429266537
    O               -4.286061502245     0.450668011111    -4.194429266537
    H                1.155774707267     1.915774341749    -3.839475176128
    H               -0.598275538017    -2.758807464732    -2.376100276707
    H               -0.598275538017    -2.758807464732     2.376100276707
    H                1.155774707267     1.915774341749     3.839475176128
    H                5.027145477405    -1.176073150085    -2.538353872066
    H                5.027145477405    -1.176073150085     2.538353872066

""")

GEOS['%s-%s-%s' % (dbse, 'P8endo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.320986373344     0.361440193807    -2.008126003438
    C               -0.253915073556     1.713282371102    -1.446792826608
    C               -0.253915073556     1.713282371102     1.446792826608
    C                2.320986373344     0.361440193807     2.008126003438
    C                2.029864844401    -2.379616120081    -1.255659867169
    C                2.029864844401    -2.379616120081     1.255659867169
    O                3.806438567173     1.351363798243     0.000000000000
    C               -2.612670939179     0.292099031874    -2.207184998419
    C               -2.612670939179     0.292099031874     2.207184998419
    N               -3.837157693467    -0.389382699166     0.000000000000
    O               -3.338960499634    -0.193375957554     4.303487249709
    O               -3.338960499634    -0.193375957554    -4.303487249709
    H                3.150360193935     0.783733056532    -3.837680523321
    H               -0.260490205070     3.591143203906    -2.292329640422
    H               -0.260490205070     3.591143203906     2.292329640422
    H                3.150360193935     0.783733056532     3.837680523321
    H                1.698591389110    -3.934948218600    -2.536959293467
    H                1.698591389110    -3.934948218600     2.536959293467
    H               -5.436813040410    -1.429496221048     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'P8exo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.601815246643     0.924566822383    -2.007405673615
    C               -0.554815123205    -1.007418463728    -1.451001293332
    C               -0.554815123205    -1.007418463728     1.451001293332
    C                1.601815246643     0.924566822383     2.007405673615
    C                4.034955492729    -0.376089962182    -1.255026098877
    C                4.034955492729    -0.376089962182     1.255026098877
    O                1.241111517692     2.673312276604     0.000000000000
    C               -3.115122822069     0.036210151668    -2.198812104876
    C               -3.115122822069     0.036210151668     2.198812104876
    N               -4.418769844038     0.555967831725     0.000000000000
    O               -3.901442524383     0.377750583007     4.300615889578
    O               -3.901442524383     0.377750583007    -4.300615889578
    H                1.497982129051     1.847378640288    -3.839612349612
    H               -0.251016959104    -2.826427461824    -2.364837705338
    H               -0.251016959104    -2.826427461824     2.364837705338
    H                1.497982129051     1.847378640288     3.839612349612
    H                5.354532696497    -1.262787533293    -2.537274357911
    H                5.354532696497    -1.262787533293     2.537274357911
    H               -6.156117945972     1.344354339031     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'P9endo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.554051073384     0.086311262521    -2.123616689415
    C               -0.973144579040     1.460936020022    -1.444377065135
    C               -0.973144579040     1.460936020022     1.444377065135
    C                1.554051073384     0.086311262521     2.123616689415
    C                1.275064912456    -2.621140678850    -1.259394974011
    C                1.275064912456    -2.621140678850     1.259394974011
    C                3.229739734700     1.139174840821     0.000000000000
    C               -3.334260245024     0.042569846594    -2.157299507688
    C               -3.334260245024     0.042569846594     2.157299507688
    O               -4.596043144515    -0.691399546664     0.000000000000
    O               -4.126162065166    -0.452121676411     4.196726175078
    O               -4.126162065166    -0.452121676411    -4.196726175078
    H                2.156326769272     0.377957142340    -4.066650862572
    H               -1.035017858088     3.336650614318    -2.298316013279
    H               -1.035017858088     3.336650614318     2.298316013279
    H                2.156326769272     0.377957142340     4.066650862572
    H                0.926351364125    -4.209802646383    -2.498541705134
    H                0.926351364125    -4.209802646383     2.498541705134
    H                3.363875453387     3.201394553009     0.000000000000
    H                5.116009212593     0.308110384532     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'P9exo', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.206896282873     0.592343149897    -2.120248754242
    C               -0.971722553203    -1.283779638545    -1.446073659060
    C               -0.971722553203    -1.283779638545     1.446073659060
    C                1.206896282873     0.592343149897     2.120248754242
    C                3.595237598010    -0.719198441627    -1.258794574595
    C                3.595237598010    -0.719198441627     1.258794574595
    C                0.902897405076     2.546589769844     0.000000000000
    C               -3.542923917860    -0.280811288707    -2.153326350494
    C               -3.542923917860    -0.280811288707     2.153326350494
    O               -4.905973268236     0.265593454406     0.000000000000
    O               -4.403918143822     0.065800017324     4.194723976040
    O               -4.403918143822     0.065800017324    -4.194723976040
    H                1.147676456511     1.257120059020    -4.064652932044
    H               -0.706926201150    -3.118492855620    -2.342746410118
    H               -0.706926201150    -3.118492855620     2.342746410118
    H                1.147676456511     1.257120059020     4.064652932044
    H                4.927680986611    -1.648192832324    -2.501396363318
    H                4.927680986611    -1.648192832324     2.501396363318
    H               -0.920807802386     3.519177058332     0.000000000000
    H                2.419882649604     3.939063378582     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'butadiene', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -3.468216365637    -0.223291847497     0.000000000000
    C               -1.146272586522     0.749237669757     0.000000000000
    C                1.146272586522    -0.749237669757     0.000000000000
    C                3.468216365637     0.223291847497     0.000000000000
    H               -3.772085595126    -2.251679349006     0.000000000000
    H               -5.132486709665     0.967856314440     0.000000000000
    H               -0.903048885942     2.790788795531     0.000000000000
    H                0.903048885942    -2.790788795531     0.000000000000
    H                3.772085595126     2.251679349006     0.000000000000
    H                5.132486709665    -0.967856314440     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'cyclohexadiene', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -2.671182948042     0.122327880670     0.443661149546
    C               -1.361349852164     0.211979091125     2.598278327326
    C                1.361349852164    -0.211979091125     2.598278327326
    C                2.671182948042    -0.122327880670     0.443661149546
    C                1.362308832887     0.467959372755    -1.999597361395
    C               -1.362308832887    -0.467959372755    -1.999597361395
    H               -4.707471738294     0.362485669716     0.446444690955
    H               -2.308773047399     0.552913024015     4.383812397059
    H                2.308773047399    -0.552913024015     4.383812397059
    H                4.707471738294    -0.362485669716     0.446444690955
    H                2.392087795467    -0.337029109670    -3.598288745833
    H                1.393432235518     2.529721193645    -2.274310457658
    H               -2.392087795467     0.337029109670    -3.598288745833
    H               -1.393432235518    -2.529721193645    -2.274310457658

""")

GEOS['%s-%s-%s' % (dbse, 'cyclopentadiene', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                0.000000000000    -2.214054512650     0.360259963025
    C                0.000000000000    -1.379363036884    -2.034813486876
    C               -0.000000000000     1.379363036884    -2.034813486876
    C               -0.000000000000     2.214054512650     0.360259963025
    C                0.000000000000     0.000000000000     2.108537296537
    H                0.000000000000    -4.161921281295     0.977530560014
    H                0.000000000000    -2.542653847531    -3.717321069475
    H               -0.000000000000     2.542653847531    -3.717321069475
    H               -0.000000000000     4.161921281295     0.977530560014
    H               -1.655029738306    -0.000000000000     3.360075385043
    H                1.655029738306     0.000000000000     3.360075385043

""")

GEOS['%s-%s-%s' % (dbse, 'ethene', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                1.249942225839     0.000000000000     0.000000000000
    C               -1.249942225839    -0.000000000000     0.000000000000
    H                2.326100952656    -1.744426289672     0.000000000000
    H                2.326100952656     1.744426289672     0.000000000000
    H               -2.326100952656    -1.744426289672     0.000000000000
    H               -2.326100952656     1.744426289672     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'ethine', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                0.000000000000     0.000000000000    -1.130604821908
    H                0.000000000000     0.000000000000    -3.142560315290
    C                0.000000000000     0.000000000000     1.130604821908
    H                0.000000000000     0.000000000000     3.142560315290

""")

GEOS['%s-%s-%s' % (dbse, 'furane', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.052191829371     0.000000000000     1.074327643117
    C                1.348667180993     0.000000000000    -1.383750141079
    C               -1.348667180993    -0.000000000000    -1.383750141079
    C               -2.052191829371    -0.000000000000     1.074327643117
    O                0.000000000000     0.000000000000     2.593814828153
    H                3.858990046905     0.000000000000     2.012181519060
    H                2.589493687518     0.000000000000    -2.999666435175
    H               -2.589493687518    -0.000000000000    -2.999666435175
    H               -3.858990046905    -0.000000000000     2.012181519060

""")

GEOS['%s-%s-%s' % (dbse, 'maleine', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.118054929941     0.000000000000     1.179813632692
    C                1.253719126326     0.000000000000    -1.486619969663
    C               -1.253719126326    -0.000000000000    -1.486619969663
    C               -2.118054929941    -0.000000000000     1.179813632692
    O                0.000000000000     0.000000000000     2.697420967246
    O                4.201094401886     0.000000000000     2.010018413293
    O               -4.201094401886    -0.000000000000     2.010018413293
    H                2.562981613715     0.000000000000    -3.051922559945
    H               -2.562981613715    -0.000000000000    -3.051922559945

""")

GEOS['%s-%s-%s' % (dbse, 'maleinimide', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                2.166949906476     0.000000000000     0.750548399187
    C                1.255230297417     0.000000000000    -1.925515598030
    C               -1.255230297417    -0.000000000000    -1.925515598030
    C               -2.166949906476    -0.000000000000     0.750548399187
    N                0.000000000000     0.000000000000     2.219080195503
    O                4.305804430778     0.000000000000     1.508177179376
    O               -4.305804430778    -0.000000000000     1.508177179376
    H                2.550797020186     0.000000000000    -3.503069239800
    H               -2.550797020186    -0.000000000000    -3.503069239800
    H                0.000000000000     0.000000000000     4.120638323033

""")

