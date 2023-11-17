module slaarnaag
use dsp
implicit none
real(dp),parameter :: pi=3.14159265358979324d0
real(dp),parameter :: twopi= 3.14159265358979324d0*2.d0
real(dp),parameter :: fourpi=3.14159265358979324d0*4.d0
real(dp),parameter :: pi_over_two=3.14159265358979324d0/2.d0
real(dp),parameter :: pi_over_four=3.14159265358979324d0/4.d0
real(dp),parameter :: threepi_over_two=1.5d0*pi
real(dp),parameter :: twopi_over_three=twopi/3.d0
real(dp),parameter :: root_pi=1.77245385090551602d0
real(dp),parameter :: root_two=1.41421356237309504d0
real(dp),parameter :: root_fourpi=3.54490770181103205546d0
real(dp),parameter :: two_root_pi=2.d0*root_pi
real(dp),parameter :: two_over_root_pi=2.d0/root_pi
real(dp),parameter :: four_over_root_pi=4.d0/root_pi
real(dp),parameter :: root_one_over_pi=0.564189583547756286d0
real(dp),parameter :: root_eight_over_pi=1.59576912160573071d0
real(dp),parameter :: root_three_over_two=0.86602540378443864676d0
real(dp),parameter :: root_threeoverfourpi=0.48860251190291992146d0
real(dp),parameter :: root_fiveoverfourpi=0.63078313050504001190d0
real(dp),parameter :: root_45_over_16pi=0.94617469575756001809d0
real(dp),parameter :: root_175_over_16pi=1.86588166295057695707d0
real(dp),parameter :: root_11025_over_256pi=3.7024941420321506331d0
real(dp),parameter :: root_43659_over_256pi=7.36787031456568657572d0
real(dp),parameter :: one_over_pi=1.d0/pi
real(dp),parameter :: one_over_twopi=1.d0/twopi
real(dp),parameter :: one_over_fourpi=1.d0/fourpi
real(dp),parameter :: one_over_root_pi=1.d0/root_pi
real(dp),parameter :: one_over_root_fourpi=1.d0/root_fourpi
real(dp),parameter :: one_over_fourpisquared=one_over_fourpi/pi
real(dp),parameter :: one_over_root_two=0.70710678118654752d0
real(dp),parameter :: half=1.d0/2.d0
real(dp),parameter :: third=1.d0/3.d0
real(dp),parameter :: twothirds=2.d0/3.d0
real(dp),parameter :: fourthirds=twothirds+twothirds
real(dp),parameter :: threefifths=3.d0/5.d0
real(dp),parameter :: sixth=1.d0/6.d0
real(dp),parameter :: sixsevenths=6.d0/7.d0
real(dp),parameter :: ninth=1.d0/9.d0
real(dp),parameter :: twoninths=2.d0*ninth
real(dp),parameter :: fourninths=4.d0*ninth
real(dp),parameter :: tenninths=1.d0+ninth
real(dp),parameter :: eleventh=1.d0/11.d0
real(dp),parameter :: twelfth=1.d0/12.d0
real(dp),parameter :: fifteenth=1.d0/15.d0
real(dp),parameter :: eighteenth=1.d0/18.d0
real(dp),parameter :: twentyoneth=1.d0/21.d0
real(dp),parameter :: fivetwentyoneths=5.d0*twentyoneth
real(dp),parameter :: thirtieth=1.d0/30.d0
real(dp),parameter :: threethirtyfifths=3.d0/35.d0
real(dp),parameter :: thirtysixth=1.d0/36.d0
real(dp),parameter :: rec_54=1.d0/54.d0
real(dp),parameter :: zero=0.d0
real(dp),parameter :: euler=0.577215664901532860606d0
real(dp),parameter :: two_euler=2.d0*euler
real(dp),parameter :: half_log_pi=0.57236494292470008706d0
real(dp),parameter :: log2=0.69314718055994530941d0
complex(dp),parameter :: czero=(0.d0,0.d0)
complex(dp),parameter :: c_one=(1.d0,0.d0)
complex(dp),parameter :: zi=(0.d0,1.d0)
logical,parameter  :: true=.true.,false=.false.
real(dp),parameter :: min_exp=1.d-150,max_exp=1.d150
real(dp) minimum_exp_arg_actual,maximum_exp_arg_actual,minimum_exp_arg&
&,maximum_exp_arg,largest_representable_number,smallest_representable_&
&number,max_rep_int
real(dp),parameter :: emach=0.111d-15,emach2=0.123d-31
real(dp) sparse_threshold
real(dp),parameter :: epsdet=1.d-150
real(dp),parameter :: epsdet_sq=epsdet*epsdet
real(dp),parameter :: speed_light_si=299792458.d0
real(dp),parameter :: elementary_charge_si=1.602176462d-19
real(dp),parameter :: proton_mass_si=1.67262158d-27
real(dp),parameter :: electron_mass_si=9.10938188d-31
real(dp),parameter :: avogadro_si=6.02214199d23
real(dp),parameter :: planck_si=6.62606876d-34
real(dp),parameter :: hbar_si=planck_si/twopi
real(dp),parameter :: amu_si=1.d-3/avogadro_si
real(dp),parameter :: mu_0_si=4.d0*pi*1.d-7
real(dp),parameter :: molar_gas_si=8.314472d0
real(dp),parameter :: boltzmann_si=molar_gas_si/avogadro_si
real(dp),parameter :: epsilon_0_si=1.d0/(mu_0_si*        speed_light_s&
&i**2)
real(dp),parameter :: fine_structure_si=                 elementary_ch&
&arge_si**2/(4.d0*pi *epsilon_0_si*hbar_si*           speed_light_si)
real(dp),parameter :: bohr              = 1.d0
real(dp),parameter :: metre             = electron_mass_si*speed_light&
&_si* fine_structure_si/hbar_si
real(dp),parameter :: centimetre        = metre*1.d-2
real(dp),parameter :: nanometre         = metre*1.d-9
real(dp),parameter :: angstrom          = metre*1.d-10
real(dp),parameter :: electron_mass     = 1.d0
real(dp),parameter :: amu               = amu_si/electron_mass_si
real(dp),parameter :: kilogram          = 1.d0/electron_mass_si
real(dp),parameter :: gram              = kilogram*1.d-3
real(dp),parameter :: aut               = 1.d0
real(dp),parameter :: second            = speed_light_si**2*fine_struc&
&ture_si **2*electron_mass_si/hbar_si
real(dp),parameter :: millisecond       = 1.d-3*second
real(dp),parameter :: microsecond       = 1.d-6*second
real(dp),parameter :: nanosecond        = 1.d-9*second
real(dp),parameter :: picosecond        = 1.d-12*second
real(dp),parameter :: femtosecond       = 1.d-15*second
real(dp),parameter :: elementary_charge = 1.d0
real(dp),parameter :: coulomb           = 1.d0/elementary_charge_si
real(dp),parameter :: hartree           = 1.d0
real(dp),parameter :: millihartree      = 1.d-3
real(dp),parameter :: electron_volt     = elementary_charge_si/ (fine_&
&structure_si**2* electron_mass_si*speed_light_si**2)
real(dp),parameter :: millielectron_volt= electron_volt*1.d-3
real(dp),parameter :: rydberg           = 0.5d0
real(dp),parameter :: millirydberg      = rydberg*1.d-3
real(dp),parameter :: joule             = 1.d0/(fine_structure_si**2* &
&electron_mass_si*speed_light_si**2)
real(dp),parameter :: erg               = joule*1e-7_dp
real(dp),parameter :: kilojoulepermole  = joule/avogadro_si*1.d3
real(dp),parameter :: kilocalpermole    = kilojoulepermole*4.184d0
real(dp),parameter :: hertz             = planck_si*joule
real(dp),parameter :: megahertz         = hertz*1.d6
real(dp),parameter :: gigahertz         = hertz*1.d9
real(dp),parameter :: terahertz         = hertz*1.d12
real(dp),parameter :: wavenumber        = hertz*speed_light_si*1.d2
real(dp),parameter :: kelvin            = boltzmann_si*joule
real(dp),parameter :: hartreebybohr     = 1.d0
real(dp),parameter :: evbyang           = electron_volt/angstrom
real(dp),parameter :: newton            = joule/metre
real(dp),parameter :: dyne              = newton*1.d-5
real(dp),parameter :: auv               = 1.0_dp
real(dp),parameter :: angperps          = angstrom/picosecond
real(dp),parameter :: angperfs          = angstrom/femtosecond
real(dp),parameter :: bohrperps         = bohr/picosecond
real(dp),parameter :: bohrperfs         = bohr/femtosecond
real(dp),parameter :: metrepersecond    = metre/second
real(dp),parameter :: hartreebybohr3    = 1.d0
real(dp),parameter :: evbyang3          = electron_volt/angstrom**3
real(dp),parameter :: pascal            = newton/metre**2
real(dp),parameter :: megapascal        = pascal*1.d6
real(dp),parameter :: gigapascal        = pascal*1.d9
real(dp),parameter :: atmosphere        = pascal*101325.027d0
real(dp),parameter :: bar               = pascal*1.d5
real(dp),parameter :: megabar           = bar*1.d6
real(dp),parameter :: invbohr           = 1.d0
real(dp),parameter :: invmetre          = 1.d0/metre
real(dp),parameter :: invnanometre      = 1.d0/nanometre
real(dp),parameter :: invangstrom       = 1.d0/angstrom
real(dp),parameter :: hartreebybohr2    = 1.d0
real(dp),parameter :: evbyang2          = electron_volt/angstrom**2
real(dp),parameter :: newtonbymetre     = newton/metre
real(dp),parameter :: dynebycentimetre  = dyne/centimetre
real(dp),parameter :: bohr3             = 1.d0
real(dp),parameter :: metre3            = metre**3
real(dp),parameter :: centimetre3       = (metre*1.d-2)**3
real(dp),parameter :: nanometre3        = (metre*1.d-9)**3
real(dp),parameter :: angstrom3         = (metre*1.d-10)**3
real(dp),parameter :: htoev=27.2113962d0
real(dp)           :: bohr_to_angstrom=0.52917715d0
character(2) periodic_table(0:92),periodic_table_nocap(0:92)
data periodic_table/ 'Wg','H ','He','Li','Be','B ','C ','N ', 'O ','F &
&','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar', 'K ','Ca','Sc','Ti','&
&V ','Cr','Mn','Fe','Co','Ni','Cu', 'Zn','Ga','Ge','As','Se','Br','Kr'&
&,'Rb','Sr','Y ','Zr', 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','S&
&n','Sb', 'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', 'Eu'&
&,'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta', 'W ','Re','Os','I&
&r','Pt','Au','Hg','Tl','Pb','Bi','Po', 'At','Rn','Fr','Ra','Ac','Th',&
&'Pa','U '/
data periodic_table_nocap/ 'wg','h ','he','li','be','b ','c ','n ', 'o&
& ','f ','ne','na','mg','al','si','p ','s ','cl','ar', 'k ','ca','sc',&
&'ti','v ','cr','mn','fe','co','ni','cu', 'zn','ga','ge','as','se','br&
&','kr','rb','sr','y ','zr', 'nb','mo','tc','ru','rh','pd','ag','cd','&
&in','sn','sb', 'te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm'&
&, 'eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta', 'w ','re','&
&os','ir','pt','au','hg','tl','pb','bi','po', 'at','rn','fr','ra','ac'&
&,'th','pa','u '/
character(12) element_names(92),element_names_caps(92)
data element_names/'hydrogen','helium','lithium','beryllium','boron','&
&carbon','nitrogen','oxygen','fluorine','neon','sodium','magnesium','a&
&luminium','silicon','phosphorus','sulphur','chlorine','argon','potass&
&ium','calcium','scandium','titanium','vanadium','chromium','manganese&
&','iron','cobalt','nickel','copper','zinc','gallium','germanium','ast&
&atine','selenium','bromine','krypton','rubidium','strontium','yttrium&
&','zirconium','niobium','molybdenum','technetium','ruthenium','rhodiu&
&m','palladium','silver','cadmium','indium','tin','antimony','telluriu&
&m','iodine','xenon','caesium','barium','lanthanum','cerium','praesody&
&mium','neodymium','promethium','samarium','europium','gadolinium','te&
&rbium','dysprosium','holmium','erbium','thulium','ytterbium','lutetiu&
&m','hafnium','tantalum','tungsten','rhenium','osmium','iridium','plat&
&inum','gold','mercury','thallium','lead','bismuth','polonium','astati&
&ne','radon','francium','radium','actinium','thorium','protactinium','&
&uranium'/
data element_names_caps /'Hydrogen','Helium','Lithium','Beryllium','Bo&
&ron','Carbon','Nitrogen','Oxygen','Fluorine','Neon','Sodium','Magnesi&
&um','Aluminium','Silicon','Phosphorus','Sulphur','Chlorine','Argon','&
&Potassium','Calcium','Scandium','Titanium','Vanadium','Chromium','Man&
&ganese','Iron','Cobalt','Nickel','Copper','Zinc','Gallium','Germanium&
&','Astatine','Selenium','Bromine','Krypton','Rubidium','Strontium','Y&
&ttrium','Zirconium','Niobium','Molybdenum','Technetium','Ruthenium','&
&Rhodium','Palladium','Silver','Cadmium','Indium','Tin','Antimony','Te&
&llurium','Iodine','Xenon','Caesium','Barium','Lanthanum','Cerium','Pr&
&aesodymium','Neodymium','Promethium','Samarium','Europium','Gadoliniu&
&m','Terbium','Dysprosium','Holmium','Erbium','Thulium','Ytterbium','L&
&utetium','Hafnium','Tantalum','Tungsten','Rhenium','Osmium','Iridium'&
&,'Platinum','Gold','Mercury','Thallium','Lead','Bismuth','Polonium','&
&Astatine','Radon','Francium','Radium','Actinium','Thorium','Protactin&
&ium','Uranium'/
integer,parameter :: number_of_elements=109
real(dp),parameter :: nuclear_mass(number_of_elements) = (/ 1.00794d0,&
&4.00260d0,6.941d0, 9.012187d0,10.811d0, 12.0107d0, 14.00674d0, 15.999&
&4d0, 18.99840d0, 20.1797d0, 22.98977d0, 24.3050d0, 26.98154d0, 28.085&
&5d0, 30.97376d0, 32.066d0, 35.4527d0, 39.948d0, 39.0983d0, 40.078d0, &
&44.95591d0, 47.867d0, 50.9415d0, 51.9961d0, 54.93805d0, 55.845d0, 58.&
&93320d0, 58.6934d0, 63.546d0, 65.39d0, 69.723d0, 72.61d0, 74.92160d0,&
& 78.96d0, 79.904d0, 83.80d0, 85.4678d0, 87.62d0, 88.90585d0, 91.224d0&
&, 92.90638d0, 95.94d0, 98.0d0, 101.07d0, 102.90550d0, 106.42d0, 107.8&
&682d0, 112.411d0, 114.818d0, 118.710d0, 121.760d0, 127.60d0, 126.9044&
&7d0, 131.29d0, 132.90545d0, 137.327d0, 138.9055d0, 140.116d0, 140.907&
&65d0, 144.24d0, 145.0d0, 150.36d0, 151.964d0, 157.25d0, 158.92534d0, &
&162.50d0, 164.93032d0, 167.26d0, 168.93421d0, 173.04d0, 174.967d0, 17&
&8.49d0, 180.9479d0, 183.84d0, 186.207d0, 190.23d0, 192.217d0, 195.078&
&d0, 196.96655d0, 200.59d0, 204.3833d0, 207.2d0, 208.98038d0, 209.0d0,&
& 210.0d0, 222.0d0, 223.0d0, 226.0d0, 227.0d0, 232.0381d0, 231.03588d0&
&, 238.0289d0, 237.0d0, 244.0d0, 243.0d0, 247.0d0, 247.0d0, 251.0d0, 2&
&52.0d0, 257.0d0, 258.0d0, 259.0d0, 262.0d0, 261.0d0, 262.0d0, 263.0d0&
&, 264.0d0, 265.0d0, 268.0d0/)
real(dp) nuclear_mass_au(number_of_elements)
real(dp),parameter :: speed_light_au=299792458.d0*metrepersecond
real(dp),parameter :: minus_1_over_2c_squared=-0.5d0/(speed_light_au**&
&2)
real(dp),parameter :: one_over_4c_squared=0.25d0/(speed_light_au**2)
real(dp),parameter :: minus_1_over_8c_squared=-0.125d0/(speed_light_au&
&**2)
end module slaarnaag
