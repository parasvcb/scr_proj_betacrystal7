# betasheetnanocrystal

The overall structure of this project and this repo will be highlighted below,

Md simulations were handled differently and after that the following process was used,

1. Data were analysed and compared for both the veloctites of buehler version and its 10 times faster velocity pull.
2. Data were traditionally analysed in terms of a) F Vs D behaviour, b) toughness values, c) A and B for the first maximal pullout force and till first peak d) hydrogen bond count at the peak of system e) the delienation of these counts being stable and stochastic and their overall frequency distribution, f) count of these system variables with the displacement.
3. Above listed two points were segregated into two different parts, data generation, and plotting, which will be detailed below.

## inspiration

The rupture of betasheet nanocrystal under the application of external force is well known. This knowledgeful wealth has largely amassed from insilico pullout experiments, where different modes of pulling from different regions and subjecting different subunits of nanocrystal has converged to highlight that perpendicularly oriented backbone H bonds in beta sheets are the major contributors to strength of silk fibers. Though not majorly, contributions of sidechain interactions in nanocrystalline strength were also acknowledged. To best of our knowledge, the nanocrystalline structure that was repeatedly used in these studies were either representatives of well conserved AGAGAG repeats or repeats of AAAAA, with exceptions of handful of structural and MD studies where doping of S,Y,W at termini were used. Relook at organization and buildup of nanocrystalline structure ignited the curiosity that even though beta sheet favoring amino acids were presnet scattered in amorphous region or region other than nanocrystalline regime. Whether this discrepancy is random or evolutionary choosen ?. More weight/ emphasis can be put on the latter as sequence repeats of them are indeed evolutionary well conserved. Hereby, in this research we would like to compare and constrast the silk nanocrystalline region with homopolymers of well known amino acids which favors and imparts stability to beta sheets in globular proteins, with exceptions of the charged and bulky amino acids. We would like to see, if further strength and stability/resilience can be imparted to these nanocrystals if alanine could have been replaced by T,I,V,N,AG,G,A.
In addition to providing the comparative mechanical basis of strength, this may also adds to helps to fine tune bio inspired materials as per initiatives of materiomics.

## Research methodology

1.  2slk structure, that does have antiparallel geometry with 5 strands organized into 3 layers were choosen (pic ref), and modifed such that each layer does have 7 strands. Original 2slk sequence repeats of 6 with two terminal modifications at N and C were modified such that, terminus modifications were removed along with mutation of glycine to alanine. This poly A structure was later mutated to other amino acids by replacing its sidechains (swapaa routine chimera) [bin/structure_preparation/7strandfrom5_logic.txt].

2.  After sicdechain modification, appropriate care was given such that clashes between the layers were minimized while having vanderwall contacts [bin/structure_preparation/resolvingclash_manually.tcl]
3.  Running MD simulations, [bin/running_mdsimulation/readme]
4.  Figure generation
    In terms of figures.

            ##### Fig 1

            is the nanocrystal representative [generated in chimera]

            ##### Fig 2

            shows interlayer sidechain packing (from transverse section) [emphasis on the central strand to be pulled (BCS)], representated using [bin/moviestuff/ftofgly.tcl]

            ##### Fig 3

            shows the profile of F vs D of one replicate, rest other are shown in supplimentary, dataframe is generated from [bin/data_generation/call3layer.py] and their cumulation is generated using [bin/plotting/vsdisp_force_fig3.R]

            ##### Fig 4

            shows mean force magnitude, mean toughness for each polymer from 3 replicates data generated from [bin/data_analysis/P_average_properties_fig4,5.py], which takes input as dir with 3 replicates and is plotted using [bin/plotting/barplotcutomized.R]. python program can be called after successfull call to [bin/data_generation/call3layer.py] for every replicate's polymers.

            ##### Fig 5

            shows mean Hb at peaks*data in variety of fashion and manners and involves follwing suit of programs
            a) [bin/plotting/barplot_hbatpeaks.R]
            This program will create two output PDFs and its sequence of running depends on the output generated after successfull call to [bin/data_generation/call3layer.py] and its cumulative analysis by the [bin/data_analysis/P_average_properties_fig4,5.py], which generates input file to this program [results/data_tsv*/out_hbonds_peak_1.tsv]
            Two output PDF's will be [Hb data at Peak1 only]
            * Fig5*segregated_plot : all, adjacent, and non adjacent categories will be in rows and tags stable, stochastic will be in columns, with color coding of hbtypes and bars of stddev and total 6 cells.
            * Fig5*raw_plot_facet: this is 3 cell (horizontal alignment), where unsegregated stable stochastic values will be used and plotted in stacked barplot (color Hb type stacking).
            b) [/bin/plotting/S_fraction_occurence_hbonds.R]
            its input files will be generated from successfull call to [bin/data_generation/call3layer.py], and there are three different types of input files were generated (kept in simulation folders and not in results section of this repo so far), (in processed folder of replicates, underscore will be preceeded by hb tags, adjacent all nonadjacent)
            * \_occurence.tsv
            <!-- hbtype identity occurence
            Main*Main SegBP1-GLY7-Main-C SegCP1-GLY2-Main-N 0.46
            Main_Main SegBP1-GLY7-Side-OT1 SegCP1-GLY2-Main-N 0.41
            Main_Main SegCP1-GLY2-Main-CA SegDP1-GLY6-Main-O 0.09
            Main_Main SegBP1-GLY5-Main-O SegCP1-GLY4-Main-N 0.48
             -->
            * \_occurence_definedbins.tsv
            <!-- hbtype Bins values FrequencyClass FrequencyTotal
            Main_Main 0_0.1 17 0.57 0.57
            Main_Main 0.1_0.2 1 0.03 0.03
            Main_Main 0.2_0.3 1 0.03 0.03
            Main_Main 0.3_0.4 2 0.07 0.07 \* \_atomdetailed.log -->
            * _atomdetailed.log

            <!--
            Hbtype totalBonds
            Main_Main 53
            Main_Side 11
            Side_Side 0

            Occ_range HB_identity hbtype
            0_0.1 SegCP1-GLY2-Main-CA : SegDP1-GLY6-Main-O M_a_i_n**_M_a_i_n
            0_0.1 SegCP1-GLY7-Main-C : SegDP1-GLY2-Main-N M_a_i_n_**M_a_i_n -->

            First will have segregation tag, along with idenity listed and overall occurenece during the frame range per HB
            Second will have the segregation tag, and occurence frequency bin range (within segregation tag and from total) along with count of Hb in that tag,
            Third one is evident.

            This program does the automatic binning of the input data of default ranges (occ of all hb's) and name them occurence (without defined bins tag), for input filetypes other than custom they have default binned segregation and this will out these pdf's
            hbond_frequency_adjacent_occurence_definedbins_subtypeFrequency.pdf
            hbond_frequency_adjacent_occurence_definedbins_wholeFrequency.pdf
            hbond_frequency_adjacent_occurence_subtypeFrequency.pdf
            hbond_frequency_all_occurence_definedbins_subtypeFrequency.pdf
            hbond_frequency_all_occurence_definedbins_wholeFrequency.pdf
            hbond_frequency_all_occurence_subtypeFrequency.pdf
            hbond_frequency_nonadj_occurence_definedbins_subtypeFrequency.pdf
            hbond_frequency_nonadj_occurence_definedbins_wholeFrequency.pdf
            hbond_frequency_nonadj_occurence_subtypeFrequency.pdf

            c) [bin/plotting/S_vsdisp_hbondcount_aspirations.R]
            above program plots hbbond dispav values with displacement and has two functions
                * plotcombined(dfgen,'hb',paste0(out,'hbcombined_vericalStack'))
                above function will only deal with the hbond categories (mc sc and combination) other than partial peak views, with their adjacent and non adjacnet and all categories
                * plotallad(dfgen,'hbadj',paste0(out,'hbAdj_verticalStack'))
                * plotallad(dfgen,'hball',paste0(out,'hbAll_verticalStack'))
                * plotallad(dfgen,'hbnad',paste0(out,'hbNad_verticalStack'))
                * above three statements will account for the biased view having stochastic stable fractions with segregations into MCSC tags but only for adjacent Hbonds (peak 1 and 2), in addition to biased view, it will also present the unbiased/unsegregated fractions and will term them as MCMC, SCSC ,MCSC for whole displacement range
            This will out PDF's
            plots/archive/plots/new4sephbAdj_verticalStackbiased1.pdf
            plots/archive/plots/new4sephbAdj_verticalStackunbiased.pdf
            plots/archive/plots/new4sephbAll_verticalStackbiased1.pdf
            plots/archive/plots/new4sephbAll_verticalStackunbiased.pdf
            plots/archive/plots/new4sephbcombined_vericalStackunbiased.pdf
            d) [bin/plotting/P_framewise_hbbonds_heatmap_tiles.R] and [bin/plotting/S_framewise_hbbonds_point.R]
            Their input files are dependent on execution of [bin/data_analysis/P_hb_integrity_cooperativity_ring_wise.py], which was manually fed the information of frames (till 1st rupture point) to be looked into and thereafter to analyse the individual Hb patterns of only 0th frame main chains HB's in form of separate Hb's, Hb rings (6) and then superrings(3). This was used to proxy the cooperativity of Hb breakage patterns.

5.  Visualization [TCL]
    TCL scripts layout (bin/moviestuff)
    TCL scripts were used to view the packing of beta strands (BCS effectively), overall layout of the nanocrsytal, rupture of H bonds with effective visualization showing only the middle stack with coplored H bond types (MCMC,MCSC,SCSC)
    following were the major type of files used:

    DCD visualization:

    1. complete_designall.tcl (/bin/plotting/images/complete_designall.png): shows colored cartoon representation with licorice of backbone (MCMC) H-bond atoms and shows effectively all Hb-types.
       [color codes: 28(magneta-MCSC), 16 (black-MCMC) , 20 (green-SCSC)]
    2. complete_designMCMC.tcl (/bin/plotting/images/complete_designMCMC.png) [(shows only MCMC in black (16)]
    3. complete_designnotMCMC.tcl (/bin/plotting/images/complete_designnotMCMC.png) [(shows both MCSC and SCSC)]
    4. designnew.tcl(/bin/plotting/images/designnew.png) and design.tcl (/bin/plotting/images/design.png) dcd and shows rupture with H bond vis and CB in vDW.

    Structure represenation:

    1. ftofgly.tcl (/bin/plotting/images/ftofgly.png) (shows the sidechain supramolecular packing of layers above and below with BCS)
    2. newdesign.tcl (/bin/plotting/images/newdesign.png) (used to show the overall nanocrystal composition and geometry of layers above and below of BCS with their sidechains)
    3. viewchains.tcl (/bin/plotting/images/viewchains.png) (serves similar purpose of newdesign.tcl but with different view)
       in menton_machine's project

6.  program running sequence
#####pending
<!--

#### pressure equilibration of systems

#### parallel run for the systems over whioch they got exlcuded was also supposed to be notified their alterations and rest

##### this segregation can be attributed to details and the readme part later.

once done with above details, here i would like to detail the data generation steps.

## Pending -->
