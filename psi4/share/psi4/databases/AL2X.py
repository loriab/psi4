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

# <<< AL2X Database Module >>>
dbse = 'AL2X'

# <<< Database Members >>>
HRXN = ['1', '2', '3', '4', '5', '6', '7', ]
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, '1'                     )] = ['%s-%s-reagent'      % (dbse, 'al2h6'),
                                                               '%s-%s-reagent'      % (dbse, 'alh3') ]
RXNM['%s-%s'            % (dbse, '1'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '1')], [-1, 2]))

ACTV['%s-%s'            % (dbse, '2'                     )] = ['%s-%s-reagent'      % (dbse, 'al2f6'),
                                                               '%s-%s-reagent'      % (dbse, 'alf3') ]
RXNM['%s-%s'            % (dbse, '2'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '2')], [-1, 2]))

ACTV['%s-%s'            % (dbse, '3'                     )] = ['%s-%s-reagent'      % (dbse, 'al2cl6'),
                                                               '%s-%s-reagent'      % (dbse, 'alcl3') ]
RXNM['%s-%s'            % (dbse, '3'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '3')], [-1, 2]))

ACTV['%s-%s'            % (dbse, '4'                     )] = ['%s-%s-reagent'      % (dbse, 'al2br6'),
                                                               '%s-%s-reagent'      % (dbse, 'albr3') ]
RXNM['%s-%s'            % (dbse, '4'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '4')], [-1, 2]))

ACTV['%s-%s'            % (dbse, '5'                     )] = ['%s-%s-reagent'      % (dbse, 'al2me4'),
                                                               '%s-%s-reagent'      % (dbse, 'alme2') ]
RXNM['%s-%s'            % (dbse, '5'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '5')], [-1, 2]))

ACTV['%s-%s'            % (dbse, '6'                     )] = ['%s-%s-reagent'      % (dbse, 'al2me5'),
                                                               '%s-%s-reagent'      % (dbse, 'alme2'),
                                                               '%s-%s-reagent'      % (dbse, 'alme3') ]
RXNM['%s-%s'            % (dbse, '6'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '6')], [-1, 1, 1]))

ACTV['%s-%s'            % (dbse, '7'                     )] = ['%s-%s-reagent'      % (dbse, 'al2me6'),
                                                               '%s-%s-reagent'      % (dbse, 'alme3') ]
RXNM['%s-%s'            % (dbse, '7'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '7')], [-1, 2]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1'                     )] =   35.8
BIND['%s-%s'            % (dbse, '2'                     )] =   52.6
BIND['%s-%s'            % (dbse, '3'                     )] =   31.4
BIND['%s-%s'            % (dbse, '4'                     )] =   28.4
BIND['%s-%s'            % (dbse, '5'                     )] =   37.8
BIND['%s-%s'            % (dbse, '6'                     )] =   30.0
BIND['%s-%s'            % (dbse, '7'                     )] =   21.4
# gmtkn24 10.1021/ct900489g Table S13
# Taken from Johnson, E. R.; Mori-Sanchez, P.; Cohen, A. J; Yang, W. J. Chem. Phys. 2008, 129, 204112.
# (back-corrected experimental values); all values are in kcal/mol.

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Binding: Al2H6 --> 2 AlH3"""
TAGL['%s-%s'            % (dbse, '2'                     )] = """Binding: Al2F6 --> 2 AlF3"""
TAGL['%s-%s'            % (dbse, '3'                     )] = """Binding: Al2Cl6 --> 2 AlCl3"""
TAGL['%s-%s'            % (dbse, '4'                     )] = """Binding: Al2Br6 --> 2 AlBr3"""
TAGL['%s-%s'            % (dbse, '5'                     )] = """Binding: Al2(CH3)4 --> 2 Al(CH3)2"""
TAGL['%s-%s'            % (dbse, '6'                     )] = """Binding: Al2(CH3)5 --> Al(CH3)2 + Al(CH3)3"""
TAGL['%s-%s'            % (dbse, '7'                     )] = """Binding: Al2(CH3)6 --> 2 Al(CH3)3"""
TAGL['%s-%s-reagent'    % (dbse, 'al2br6'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'al2cl6'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'al2f6'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'al2h6'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'al2me4'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'al2me5'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'al2me6'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'albr3'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alcl3'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alf3'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alh3'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alme2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'alme3'                 )] = """ """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'al2br6', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     3.174154631687
    AL               0.000000000000     0.000000000000    -3.174154631687
    BR              -0.000000000000     3.355131304663     0.000000000000
    BR               0.000000000000    -3.355131304663     0.000000000000
    BR              -3.705165443212    -0.000000000000     5.251115673995
    BR               3.705165443212     0.000000000000     5.251115673995
    BR               3.705165443212     0.000000000000    -5.251115673995
    BR              -3.705165443212    -0.000000000000    -5.251115673995

""")

GEOS['%s-%s-%s' % (dbse, 'al2cl6', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     3.014941633736
    AL               0.000000000000     0.000000000000    -3.014941633736
    CL              -0.000000000000     3.067175307365     0.000000000000
    CL               0.000000000000    -3.067175307365     0.000000000000
    CL              -3.441622534928    -0.000000000000     4.921401693216
    CL               3.441622534928     0.000000000000     4.921401693216
    CL               3.441622534928     0.000000000000    -4.921401693216
    CL              -3.441622534928    -0.000000000000    -4.921401693216

""")

GEOS['%s-%s-%s' % (dbse, 'al2f6', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     2.602014367842
    AL               0.000000000000     0.000000000000    -2.602014367842
    F               -0.000000000000     2.207272302695     0.000000000000
    F                0.000000000000    -2.207272302695     0.000000000000
    F               -2.724072880069    -0.000000000000     4.055439357271
    F                2.724072880069     0.000000000000     4.055439357271
    F                2.724072880069     0.000000000000    -4.055439357271
    F               -2.724072880069    -0.000000000000    -4.055439357271

""")

GEOS['%s-%s-%s' % (dbse, 'al2h6', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     2.469047470458
    AL               0.000000000000     0.000000000000    -2.469047470458
    H               -0.000000000000     2.164667217662     0.000000000000
    H                0.000000000000    -2.164667217662     0.000000000000
    H               -2.665875895164    -0.000000000000     3.773354010748
    H                2.665875895164     0.000000000000     3.773354010748
    H                2.665875895164     0.000000000000    -3.773354010748
    H               -2.665875895164    -0.000000000000    -3.773354010748

""")

GEOS['%s-%s-%s' % (dbse, 'al2me4', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                3.323139154250     0.000000000000     4.142660340781
    AL               0.000000000000     0.000000000000     2.482575906157
    C               -3.323139154250    -0.000000000000     4.142660340781
    H               -0.000000000000     2.169643190952     0.000000000000
    AL               0.000000000000     0.000000000000    -2.482575906157
    C                3.323139154250     0.000000000000    -4.142660340781
    H                0.000000000000    -2.169643190952     0.000000000000
    C               -3.323139154250    -0.000000000000    -4.142660340781
    H               -3.553589040873    -1.658802445830    -5.348913542081
    H               -4.874857000724    -0.000000000000    -2.785150567745
    H               -3.553589040873     1.658802445830    -5.348913542081
    H                3.553589040873    -1.658802445830    -5.348913542081
    H                3.553589040873     1.658802445830    -5.348913542081
    H                4.874857000724     0.000000000000    -2.785150567745
    H               -3.553589040873    -1.658802445830     5.348913542081
    H               -3.553589040873     1.658802445830     5.348913542081
    H               -4.874857000724    -0.000000000000     2.785150567745
    H                3.553589040873    -1.658802445830     5.348913542081
    H                4.874857000724     0.000000000000     2.785150567745
    H                3.553589040873     1.658802445830     5.348913542081

""")

GEOS['%s-%s-%s' % (dbse, 'al2me5', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                0.696514781275     3.216547172040     4.306423114286
    AL               0.460942715956    -0.011923864917     2.478128252860
    C                0.511223409627    -3.385964995521     4.044542679217
    C               -2.758376427663     0.369291173940     0.000000000000
    AL               0.460942715956    -0.011923864917    -2.478128252860
    C                0.696514781275     3.216547172040    -4.306423114286
    H                2.593031066584    -0.018694656528     0.000000000000
    C                0.511223409627    -3.385964995521    -4.044542679217
    H                2.299110586693    -3.742492925831    -5.012633800285
    H                0.252029071779    -4.896286943527    -2.662879967607
    H               -0.982166700511    -3.582169481794    -5.458198261742
    H                2.406693609881     3.312348537245    -5.457922549495
    H               -0.904057547926     3.468424711925    -5.587306125427
    H                0.721933192431     4.835243857592    -3.027796660171
    H                2.299110586693    -3.742492925831     5.012633800285
    H               -0.982166700511    -3.582169481794     5.458198261742
    H                0.252029071779    -4.896286943527     2.662879967607
    H                2.406693609881     3.312348537245     5.457922549495
    H                0.721933192431     4.835243857592     3.027796660171
    H               -0.904057547926     3.468424711925     5.587306125427
    H               -3.701289675481    -0.558978891385     1.597619570048
    H               -3.356521526368     2.339909130936     0.000000000000
    H               -3.701289675481    -0.558978891385    -1.597619570048

""")

GEOS['%s-%s-%s' % (dbse, 'al2me6', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -0.000000000000     4.558649222805     3.062953272350
    AL              -0.000000000000     2.470940756639    -0.017844055316
    C               -0.000000000000     3.890048773000    -3.464310482908
    C                3.197999115274     0.000000000000     0.421400568105
    AL               0.000000000000    -2.470940756639    -0.017844055316
    C                0.000000000000    -4.558649222805     3.062953272350
    C               -3.197999115274    -0.000000000000     0.421400568105
    C                0.000000000000    -3.890048773000    -3.464310482908
    H               -1.656342947220    -5.080491018962    -3.792843083207
    H                0.000000000000    -2.434132434563    -4.927100327665
    H                1.656342947220    -5.080491018962    -3.792843083207
    H               -1.657610770891    -5.787951132159     3.145400491870
    H                1.657610770891    -5.787951132159     3.145400491870
    H                0.000000000000    -3.415608870886     4.781815924881
    H               -1.656342947220     5.080491018962    -3.792843083207
    H                1.656342947220     5.080491018962    -3.792843083207
    H               -0.000000000000     2.434132434563    -4.927100327665
    H               -1.657610770891     5.787951132159     3.145400491870
    H               -0.000000000000     3.415608870886     4.781815924881
    H                1.657610770891     5.787951132159     3.145400491870
    H               -3.751092538223    -0.000000000000     2.405583099670
    H               -4.153434504501    -1.602157521269    -0.483806408221
    H               -4.153434504501     1.602157521269    -0.483806408221
    H                3.751092538223     0.000000000000     2.405583099670
    H                4.153434504501     1.602157521269    -0.483806408221
    H                4.153434504501    -1.602157521269    -0.483806408221

""")

GEOS['%s-%s-%s' % (dbse, 'albr3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     0.000000000000
    BR               2.114427655041     3.662296127459     0.000000000000
    BR               2.114427655041    -3.662296127459     0.000000000000
    BR              -4.228855310081     0.000000000000     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'alcl3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     0.000000000000
    CL               1.961729889657     3.397815839612     0.000000000000
    CL               1.961729889657    -3.397815839612     0.000000000000
    CL              -3.923459779313     0.000000000000     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'alf3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     0.000000000000
    F                1.545001848671     2.676021699686     0.000000000000
    F                1.545001848671    -2.676021699686     0.000000000000
    F               -3.090003697342     0.000000000000     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'alh3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    AL               0.000000000000     0.000000000000     0.000000000000
    H                1.491483526275     2.583325246161     0.000000000000
    H                1.491483526275    -2.583325246161     0.000000000000
    H               -2.982967052551     0.000000000000     0.000000000000

""")

GEOS['%s-%s-%s' % (dbse, 'alme2', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C                0.000000000000    -3.269241990451     0.458463952380
    AL               0.000000000000     0.000000000000    -1.317874540410
    C               -0.000000000000     3.269241990451     0.458463952380
    H                0.000000000000     0.000000000000    -4.318699207834
    H                0.000000000000    -3.087740551782     2.510393098953
    H                1.650567401518    -4.391277498808    -0.075285088605
    H               -1.650567401518    -4.391277498808    -0.075285088605
    H               -1.650567401518     4.391277498808    -0.075285088605
    H                1.650567401518     4.391277498808    -0.075285088605
    H               -0.000000000000     3.087740551782     2.510393098953

""")

GEOS['%s-%s-%s' % (dbse, 'alme3', 'reagent')] = qcdb.Molecule("""
    units Bohr
    no_com
    no_reorient
    0 1
    C               -0.084786254116     3.727451483689     0.000000000000
    AL               0.000000000000     0.000000000000     0.000000000000
    C                3.270460803306    -1.790298691888     0.000000000000
    C               -3.185674549190    -1.937152791801     0.000000000000
    H                1.785509789000     4.592499159707     0.000000000000
    H               -1.103964191112     4.439395481205     1.650739698801
    H               -1.103964191112     4.439395481205    -1.650739698801
    H                4.396611359726    -1.263636706231    -1.650739698801
    H                4.396611359726    -1.263636706231     1.650739698801
    H                3.084466044665    -3.842546415834     0.000000000000
    H               -4.869975833665    -0.749952743874     0.000000000000
    H               -3.292647168614    -3.175758774974     1.650739698801
    H               -3.292647168614    -3.175758774974    -1.650739698801

""")

