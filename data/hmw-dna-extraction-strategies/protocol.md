---
layout: page
title: HMW DNA extraction strategies for complex metagenomes
modified: 2021-10-21
excerpt: "by Trigodet and Lolans et al 2021"
comments: true
authors: [karen]
image:
  feature: https://merenlab.org/data/hmw-dna-extraction-strategies/images/header.png
  display: true
---

<div class="extra-info" markdown="1">

<span class="extra-info-header">Summary</span>

The intention of this page is to offer detailed descriptions of the high molecular weight DNA extraction protocols explained in our [paper by Trigodet and Lolans et al](https://doi.org/10.1111/1755-0998.13588). The sections in this document will detail the four methods covered in the paper, namely, GT (or M4), PC (or M5), and AE (or M6).

You can also find our reproducible data products and the bioinformatics workflow that underlie our key findings in the study [here](/data/hmw-dna-extraction-strategies/).

</div>

<div class="pub_float">
<div class="altmetric-embed" data-badge-type="donut" data-doi="10.1111/1755-0998.13588"></div>
<div class="__dimensions_badge_embed__" data-doi="10.1111/1755-0998.13588" data-hide-zero-citations="true" data-legend="hover-bottom" data-style="small_circle"></div>
    <span class="pub-title"><a href="https://doi.org/10.1111/1755-0998.13588" target="_new">High molecular weight DNA extraction strategies for long-read sequencing of complex metagenomes</a></span>
    <span class="pub-authors"><span class="pub-member-author">Trigodet F</span>, <span class="pub-member-author">Lolans K</span>, Fogarty E, <span class="pub-member-author">Shaiber A</span>, Morrison HG, Barreiro L, Jabri B, <span class="pub-member-author">Eren AM</span></span>
    <div class="pub-info">
    <div class="pub-featured-image">
    <a href="/images/pubs/Trigodet_and_Lolans_et_al_HMW.png"><img src="/images/pubs/Trigodet_and_Lolans_et_al_HMW.png" style="max-width: 100px; max-height: 80px; width: auto; border: none; height: auto; margin: 0 auto; display: block; transform: translateY(15%);" /></a>
    </div>
    <div class="pub-highlights">
    <span style="display: inline-block; padding-bottom: 5px;">- A study that <b>benchmarks six high molecular weight DNA extraction strategies</b> (commercially available kits, phenol-chloroform extraction, and agarose encasement followed by agarase digestion) <b>for long-read sequencing of metagenomes</b> with MinION.</span><br /><span style="display: inline-block; padding-bottom: 5px;">- It turns out the protocol that works best for sequencing DNA from microbial isolates may not be the most effetive method for long-read sequencing of metagenomes ¯\_(ツ)_/¯</span>
    </div>
    </div>
    <span class="pub-journal"> 📚 <b>Molecular Ecology Resources</b>, 22(5):1786-1802 | 🔍 <a href="http://scholar.google.com/scholar?hl=en&amp;q=High+molecular+weight+DNA+extraction+strategies+for+long-read+sequencing+of+complex+metagenomes" target="_blank">Google Scholar</a> | 🔗 <a href="https://doi.org/10.1111/1755-0998.13588" target="_blank">doi:10.1111/1755-0998.13588</a></span>
</div>


## Qiagen Genomic Tip 20/G, with enzymatic treatment,  Extraction Protocol [aka, GT or M4]

<div class="extra-info" markdown="1">

<span class="extra-info-header">PURPOSE</span>

The intention of this protocol is to isolate high molecular weight DNA from tongue dorsum samples.

Best practies include avoiding pipetting without use of wide-bore or cut off pipette tips and limiting vortexing, using mixer shakers or anything else which may generate a velocity gradient that could shear the DNA. In addition, care should be taken to limit the introduction of nucleases when making reagents, buffers, etc  through the sole use of  nuclease-free components. Avoid unnecessary heating and do not freeze isolated DNA.

Isolated DNA should be stored in the fridge. A good extraction will be stable for months.

</div>

### Reagents

#### Purchased reagents

1. RNase, Sigma-Aldrich, Cat#: R6513 – 10mg [store @ -20°C]
2. Proteinase K, Sigma-Aldrich, Cat#: P-2308-10mg [store @ -20°C]
3. Lysozyme, Sigma-Aldrich, Cat#: L6876- 1g [store @ -20°C]
4. Lysostaphin, Sigma-Aldrich, Cat#: SAE0091 - 2mg [store @ -20°C]
5. Mutanolysin, Sigma-Aldrich, Cat#: SAE0092 - 10KU [store @ -20°C]
6. Genomic DNA Buffer Set, Qiagen, Cat#: 19060
7. Genomic -tip 20/G, Qiagen, Cat#: 10223

#### Prepared Reagents

1. **RNase**: CAS: 9001-99-4. To prepare a 10-mg/mL stock, combine the following, then vortex, aliquot to 250-ul per tube, and store aliquots @ -20C:

    - 10-mg lyophilized powder
    - 75-ul 5M NaCl sterile stock
    - 100-ul 0.5M Tris-HCl, pH 7.6 sterile stock
    - add nuclease-free water to total volume of 1-ml


2. **Proteinase K**: CAS: 39450-01-6. To prepare a 20-mg/mL stock, combine the following, then vortex, and store aliquots @ -20C:

    - 10-mg lyophilized powder
    - 500-ul nuclease-free H20

3. **Mutanolysin**. To prepare a 10KU/mL stock, combine the following, then vortex, and store aliquots @ -20C::

    - 10KU lyophilized powder
    - 1000-ul molecular grade TE buffer

4. **Lysozyme**: CAS 12850-88-3. To prepare a 100mg/mL stock, combine the following, then vortex, aliquot to 250-ul per tube, and store aliquots @ -20C:

    - 100-mg powder
    - 1000-ul molecular grade TE buffer

5. **Lysostaphin**. To prepare an 8mg/mL (or 4000U/mL) stock, combine the following in the lysostaphin bottle, then vortex, and transfer to a single microfuge tube & store @ -20C:

    - 2-mg powder (= entire bottle)
    - 250-ul nuclease-free H20


### Protocol

#### Sample-type Specific Initial Steps for Tongue Dorsum Samples

{:.notice}
The purpose of the steps below relates to the fact that Tongue Dorsum samples were collected into a volume (200-500ul) of phosphate buffered saline (PBS). If pelleting of bacteria by centrifugation is not performed, then the final enzyme concentrations in subsequent steps will be askew!!

a. Vortex sample until uniform suspension

b. Centrifuge, 10,000 RPM for 10min to pellet microbial cells to bottom of tube (tabletop microfuge is sufficient)

c. Carefully decant/remove supernatant from each of the tubes. Save pellet.

d. Resuspend the pellet in 1-mL `Buffer B1`, into which 20-ul of `RNase (10mg/mL)` was added, and *vortex vigorously*.

{:.warning}
The purpose of vortexing here is to thoroughly resuspend the cells to achieve a homogenous DNA solution; DNA will not be damaged as it is still packaged inside of the cell at this point

| Step when Added | Lytic Enzyme | Stock Conc. | Enzyme stock vol. to add to 1-mL Buffer B1 | Component Added? |
| -------- | -------- | -------- | -------- | -------- |
| sample-specific protocol    | RNaseA    | 10-mg/mL     |20-uL     |     |

Continue to general procedure.

#### General Procedure

This procedure is based on the Qiagen protocol, 'Preparation of Gram-Negative and some Gram-Positive Bacterial Samples', which has been modified to include lysostaphin and mutanolysin

1. To each sample, add the additional lytic enzymes that are outlined in the table below, then mix by inversion (6-8X), and **do not vortex**:

    | Step when Added | Lytic Enzyme | Stock Conc. | Enzyme stock vol. to add to 1-mL Buffer B1 |
    | -------- | -------- | -------- | -------- |
    | general protocol, step1    | ProteinaseK    | 20-mg/mL     |45-uL     |
    general protocol, step1     |Lysozyme    | 100-mg/mL     |20-uL     |
    general protocol, step1     | Lysostaphin   | 8-ug/mL (aka 4000 U/mL)     |9-uL     |
    general protocol, step1     | Mutanolysin    | 10-KU/mL     |45-uL     |

2. Incubate @ **37°C** for **2.5 hrs** (used 37C incubator)

3. Add **0.35-mL** of **Buffer B2**, then mix by inversion (end-over-end rotation 10X), and **DO NOT VORTEX!!!!**.

    {:.notice}
    Note: It is important for the solution to be mixed thoroughly at this stage – but  do it gently.

4. Incubate @ **50°C** for **1.5 hrs** (used water bath). Some notes regarding this step:

    * We are using an extended incubation to ensure that the lysate clears.

    * The lysate **must** clear.

    * If particulate matters remains, pellet by centrifugation [**10 min at 5,000 x g**] – otherwise this matter will clog the filter in subsequent steps.

    * You don’t want to be here for 6hrs waiting for solutions to empty by gravity flow. But centrifugation, while it may be necessary, is not optimal for HMW DNA recovery.

5. Place **Buffer QF** in the **50C** water bath to warm. This will increase yields during the elution step.

6. Equilibrate a Qiagen Genomic Tip 20/G with **1-mL** of **Buffer QBT**. Allow to empty by gravity flow.

7. Once incubation step (#4 above) is complete, mix the samples of interest by **inversion** (10X end-over-end rotation). Transfer the sample to the equilibrated Genomic tip. Use **WIDE BORE PIPET TIP**, and allow to enter the resin by gravity flow.

8. Wash the Genomic tip with **3 rounds** of **1-mL Buffer QC** [i.e.   3 x 1-mL]. Change tips between Buffer QC addition between different samples, and allow Buffer QC to move through the tip by gravity flow.

9. Elute the Genomic DNA with **2 rounds** of **1-mL of warmed Buffer QF** [i.e.   2 x 1-mL], Elute into a 14-mL Falcon tube. Cannot use microfuge tubes due to downstream volumes added.

10. Precipitate the DNA by adding **1.4-mL (0.7 volumes)** room-temperature **isopropanol** to the eluted DNA. Invert the tube **10-20X** (end-over-end rotation) to GENTLY mix..

11. Recover the Precipitated DNA.

    Is there visible mass of DNA?

    1. Spool the DNA with a glass rod (a 1-mL pipet tip works pretty well too) -- once the DNA is hooked onto the rod or pipet tip, very slowly (i.e. snail pace) remove the DNA from the tube, thus, allowing the isopropanol to drip off.

    2. Transfer the spooled DNA to a **low bind microfuge tube** containing **100-ul buffer** (Buffer = 10 mM Tris-HCl pH 8-8.5).

    3. Incubate without mixing at **4°C** for **24-48hrs** (TIME:_______) to allow the pellet to fully resuspend into a translucent viscous gel.

    If there is no visible mass of DNA?

    1. Centrifuge immediately **5000 x g, 15 min**.

    2. Carefully remove the supernatant (using a pipet) without disturbing the pelleted DNA.

    3. Wash the DNA with **1-mL of cold 70% Ethanol** (Use freshly made ethanol (i.e. made day of extraction)).

    4. Mix gently (**DO NOT VORTEX**)

    5. Centrifuge immediately **5000 x g, 10min**.

    6. Carefully remove the supernatant (using a pipet) without disturbing the pelleted DNA.

    7. Air-dry for **5-10 min**

    8. Resuspend DNA in **100ul** buffer (buffer = 10 mM Tris-HCl pH 8-8.5)

    9. Transfer to a **low-bind microfuge tube**

    10. Incubate without mixing at **4°C** for **24-48hrs** (TIME:_______) to allow the pellet to fully resuspend into a translucent viscous gel.

<div class="extra-info" markdown="1">
<span class="extra-info-header">Final notes</span>
* When attempting to pipet the DNA for any subsequent downstream workflows--- you will need to use wide-bore pipet tips and pipet very slowly. It will not be unusual for the DNA to pull out of the pipet tip as it is very viscous & stringy and hopefully, long!!

* This protocol has been used with other low biomass sample types, like mucosal brushes, with slight modifications.
</div>


## Phenol/Chloroform HMW DNA Extraction Procedure [aka, PC or M5]

{:.notice}
The protocol below is adapted from a protocol outlined [here](https://www.protocols.io/view/ultra-long-read-sequencing-protocol-for-rad004-mrxc57n) in Molecular Cloning by Sambrook and Russell, 3rd edition.

<div class="extra-info" markdown="1">

<span class="extra-info-header">PURPOSE</span>

The intention of this protocol is to isolate high molecular weight DNA.

Best practies include avoiding pipetting without use of wide-bore or cut off pipette tips and limiting vortexing, using mixer shakers or anything else which may generate a velocity gradient that could shear the DNA. In addition, care should be taken to limit the introduction of nucleases when making reagents, buffers, etc  through the sole use of  nuclease-free components. Avoid unnecessary heating and do not freeze isolated DNA.

Isolated DNA should be stored in the fridge. A good extraction will be stable for months.

</div>

### Reagents

#### Purchased Reagents

1. RNase, Sigma-Aldrich, Cat#: R6513 – 10mg [store @ -20°C]
2. Proteinase K, Sigma-Aldrich, Cat#: P-2308-10mg  [store @ in -20C]
3. BioUltra TE-saturated phenol, Sigma Aldrich, Cat#: 77607 [light-sensitive, store @ 4C]
4. BioUltra Chloroform-Isoamyl Alcohol (24:1), Sigma Aldrich, Cat#: 25666
5. Ethanol Solution 96%, Molecular Biology Grade, Fisher Scientific, Cat#: BP8202500
6. 5M Ammonium Acetate (100mL), Fisher Scientific, Cat#: 50-103-5191 [note manufacturer = GBiosciences]
7. 10X PBS, pH 7.4, Fisher-Scientific, cat#: 70-011-044
8. Phase-lock gel tubes [light], QuantaBio, VWR, cat#: 10847-800
9. Ultrapure 0.5M EDTA, pH8.0, Fisher-Scientific, cat#: 15575020

#### Prepared Reagents

* **RNase**: CAS: 9001-99-4. To prepare a 10-mg/mL stock, combine the following, then vortex, aliquot to 250-ul per tube, and store aliquots @ -20C:

    - 10-mg lyophilized powder
    - 75-ul 5M NaCl sterile stock
    - 100-ul 0.5M Tris-HCl, pH 7.6 sterile stock
    - add nuclease-free water to total volume of 1-ml

*  **Proteinase K**: CAS: 39450-01-6. To prepare a 20-mg/mL stock, combine the following, then vortex, and store aliquots @ -20C:

    - 10mg lyophilized powder
    - 500ul nuclease-free H20

* **1X PBS**

* **TLB** (Total Lysis Buffer), combine the following, and add RNaseA fresh (i.e., immediately before use):

    - 100mM NaCl
    - 10mM Tris-HCL, pH 8.0
    - 25mM EDTA, pH 8.0
    - 0.5% (w/v) SDS
    - RNase A (add 20-ul of 10-mg/mL stock for every 10-mL TLB needed)

* **EB** (Elution Buffer). Solution EB (10mM Tris-HCL, pH8.5, Qiagen) can be used here [or other similarly formulated solution].

### Sample-type Specific Initial Steps for Tongue Dorsum Samples

1. Starting material used: 500-ul tongue dorsum sample
2. Use alcohol-resistant marker when labelling Falcon tubes

### Procedure

Before starting, please double-check the following notes:

* If phase-lock gel is only available in 2-ml tubes, you will need to transfer it into 15-mL conical tubes -- will need four (4) conical tubes per sample.
* Transfer can be performed by cutting the lid off one 2-ml tube, stuff if lightly into a 15-ml conical and then centrifuge it [3000xg, 3 min]. Repeat so that three (3) 2-mL tubes are combined into each conical tube. This is best done PRIOR to needing these.
* Spin green microfuge tubes 1st in microfuge to consolidate gel at tube bottom (2 min at max speed).
* Cut phase lock tube at the hinge (close to lid part) so that a lip remains. This ensures the green tubes wedge into the conical tube and that they then don’t slide inside of the conical tube when further centrifuged.

Here are the procedural steps, followed by notes, tips, or warnings:

1. To each sample, add 10-ml TLB and vortex at full speed for 5 seconds.

    - Ensure you have sufficient volume of TLB available (i.e. 10-ml TLB times the number of samples to be processed)
    - Make sure to add RNAse A to TLB before use. Add 20-ul of 10-mg/mL stock added to every 10-mL TLB needed
    - The purpose of vortexing here is to thoroughly resuspend the cells to achieve a homogenous DNA solution; DNA will not be damaged as it is still packaged inside of the cell at this point.

2. Incubate at 37°C for 1 hour. Solution will become transparent as the cells lyse.

3. Add ProteinaseK (stock solution: 20-mg/mL) to a final concentration of 200-µg/ml.

    - Add 100-ul of ProtK stock solution added. check math if any volume changes
    - Mix by slowly rotating end-over-end 3 times.

4. Incubate at 50°C for 2 hours

    - Mix every 30 minutes by slowly rotating end-over-end 3 times.
    - sample-specific note: at end of incubation, the tongue dorsum still showed particulate matter.

5. Add light phase-lock gel to 2 x 15 ml Falcons.

    - If phase-lock gel is only available in 2 ml tubes, transfer it by cutting the lid off 3 x 2 ml tubes and spinning it out into each 15 ml Falcon.
    - 15 ml conicals are used as they are narrower -- thereby decreasing the surface area of the interface/gel. Using two tubes per sample means they balance each other in the centrifuge while also providing the phenol sufficent space to move which ultimately improves the emulsion.

6. Split the viscous lysate (from step#4) into the two 15 ml phase-lock gel Falcon tubes (step#5).

    - his is easiest using a 10 ml serological pipette and pipetting at slow speed.
    - GO VERY SLOWLY -- LYSATE IS VERY VISCOUS

7. Add 5-ml recently opened BioUltra TE-saturated phenol to each Falcon tube containing lysate.

    - PERFORM THIS STEP IN CHEMICAL HOOD.

8. Place on a HulaMixer (or other gentle sample mixer) at 20 rpm for 10 minutes.

    - Use 20 rpm for end-over-end rotation
    - Note: if a fine emulsion has not formed after a minute, gradually increase the rotation speed.

9. Spin in a centrifuge at 3600 x g for 10 minutes.

10. Pour the aqueous phases into two new 15 ml Falcon tubes containing phase-lock gel.

    - Try to avoid transferring any protein which may form a white layer above the phase-lock gel.
    - Note: It helps to have plenty of light when pouring off the aqueous phase; be aware that phenol may break through -- will need to re-centrifuge if this happens.

11. Add 2.5-ml buffer saturated phenol and 2.5-ml chloroform-isoamyl alcohol 24:1 to each tube.

    - MAKE SURE CAPS ARE TIGHT!!!!!

12. Place on a HulaMixer (or other gentle sample mixer) at 20 rpm for 10 minutes.

    - If a fine emulsion has not formed after a minute gradually increase the rotation speed.
    - Make sure solutions of 100% ethanol are on ice as you will need in step#16 & beyond.

13. Centrifuge at (3600 x g) for 10 minutes.

14. Combine the aqueous phases from the two tubes by pouring slowly into a new 50 ml Falcon tube.

    - Work quickly at pouring as the phase lock will likely move.

15. Add 4-ml 5 M ammonium acetate.

16. Add 30 ml ice-cold ethanol (absolute ethanol, 100%) and watch the DNA precipitate.

    - Depending on amount of DNA -- you could see small bubbles that will, over time, pull the mass of DNA to the surface. This may appear to look like a small jellyfish with tentacles hanging down.

17. Incubate at -20C for 30-min (or longer) before next step.

    - **Do not** incubate at -80C as the extremely low temperature precipitates the salts and freezes the sample.

18. Make 70% EtOH [1-mL needed per sample]. Store 70% EtOH at -20C.

19. After -20C incubation (step#17), centrifuge the 50-mL conical at 3600 x g for 10min.

20. Carefully pour off the supernatant

21. GENTLY resuspend pelleted DNA in 1 ml ice-cold 70% ethanol & transfer to a microfuge tube. Transfer to microfuge with WIDE-BORE pipet tips.

22. Spin down at 10,000 x g (10,600 RPM x 2min) in microcentrifuge,  then carefully remove as much of the ethanol as possible.

    - Time is per Sambrook & Russell

23. Wash again with 1 ml 70% ethanol.

24. Spin down at 10,000 x g (10,600 RPM x 2min) in microcentrifuge, then remove as much of the ethanol as possible.

25. Let the remaining ethanol evaporate by leaving open at room temperature (~22-25C) for 5-10 min.

26. Add 85μl Buffer EB (Qiagen) and incubate without mixing at 4°C for 24-48hrs (TIME:_______) to allow the pellet to fully resuspend into a translucent viscous gel.

    - Modified from protocols.io version to NOT include detergents or surfactants as outlined (i.e. TritonX-100) as the nanopore protocol specifically stated to exclude these.
    - Warning: From this step forward, in subsequent downstream workflows, all pipetting steps must utilize wide bore or cut-off pipet tips.

## Agarose Plug Encasement and Extraction HMW DNA Procedure [aka, AE or M6]

{:.notice}
The protocol below is adapted from Pulsed-Field Gel Electrophoresis (PFGE) based DNA extraction Protocol of Matushek et al. 1996. Journal of Clinical Microbiology.

<div class="extra-info" markdown="1">

<span class="extra-info-header">PURPOSE</span>

The intention of this protocol is to isolate high molecular weight (HMW) DNA using the agarose plug methodology first developed for DNA assessment via PFGE.

Best practies include avoiding pipetting without use of wide-bore or cut off pipette tips and limiting vortexing, using mixer shakers or anything else which may generate a velocity gradient that could shear the DNA. In addition, care should be taken to limit the introduction of nucleases when making reagents, buffers, etc  through the sole use of  nuclease-free components. Avoid unnecessary heating and do not freeze isolated DNA.

Isolated DNA should be stored in the fridge. A good extraction will be stable for months.

</div>

### Reagents

#### Purchased Reagents

1. Sodium Chloride
2. Trizma pre-set crystals, Sigma-Aldrich, cat#: T9743-100G
3. 0.5M EDTA, pH8.0, Fisher Scientific, cat#: 15575020
4. Sodium Deoxycholate, Sigma-Aldrich, cat#: D6750-25G
5. Brij 58, Sigma-Aldrich, cat#: P5884-100G
6. N-Lauroylsarcosine sodium salt, Sigma-Aldrich, cat#: L5125-50G
7. Thermo Scientific™ TopVision Low Melting Point Agarose, Fisher Scientific, cat#: FERR0801
8. RNase, Sigma-Aldrich, Cat#: R6513 – 10mg [store @ -20°C]
9. Proteinase K, Sigma-Aldrich, Cat#: P-2308-10mg  [store @ in -20C]
10. Lysozyme, Sigma-Aldrich, Cat#: L6876- 1g [store @ -20°C]
11. Lysostaphin, Sigma-Aldrich, Cat#: SAE0091 - 2mg [store @ -20°C]
12. Mutanolysin, Sigma-Aldrich, Cat#: SAE0092 - 10KU [store @ -20°C]
13. β-Agarase I, Fisher Scientific, cat#: 50-811-726 [store @ -20°C]

#### Reagent Preparation

Following sections describe chemicals for three classes.

##### (1) Cell Lysis Chemicals

1. 5M NaCl (CAS: 7647-14-15):

    - 146-gm NaCl
    - 400-ml MilliQ H20

    Stir until dissolved. q.s. with water to 500-ml. Autoclave & store at room temp (RT)

2. 0.5M TRIS:

    - 37.25-gm TRIZMA preset crystals, pH 7.6
    - 400-ml MilliQ H20

    Stir until dissolved. q.s. with water to 500ml. Autoclave & store at RT.

3. 0.5M EDTA, pH8.0: purchased commercially

4. 10% Deoxycholate (CAS:302-94-4):

    - 25-gm deoxycholate
    - 250-ml water

    Filter sterilize (0.22uM). Store at RT.

5. 10% Brij 58 (CAS:9004-95-9):

    - 25-gm Brij 58
    - 250-ml MilliQ water

    Stir overnight (note: stock is cloudy). Do not sterilize. Store at RT.

6. 20% SLS (Sodium-N-lauroyl sarcosine) (CAS: 137-16-6):

    - 50-gm SLS
    - 250-ml MilliQ water

    Stir overnight. Filter sterilize (0.22uM). Store at RT

7. 10% SDS (Sodium lauryl sulfate) (CAS: 151-21-3):

    - 25-gm SDS
    - 250-ml water

    Stir overnight. Filter sterilize (0.22uM). Store at RT.

##### (2) Cell Lysis Enzymes

1. RNase: CAS: 9001-99-4

    To prepare a 10-mg/mL stock, combine the following:

    - 10-mg lyophilized powder
    - 75-ul 5M NaCl sterile stock
    - 100-ul 0.5M Tris-HCl, pH 7.6 sterile stock
    - add nuclease-free water to total volume of 1-ml

    Vortex. Then aliquot to 250-ul per tube. Store aliquots @ -20C

2. Proteinase K: CAS: 39450-01-6

    To prepare a 20-mg/mL stock, combine the following:

    - 10-mg lyophilized powder
    - 500-ul nuclease-free H20

    Vortex. Store aliquots @ -20C

3. Mutanolysin:

    To prepare a 10KU/mL stock, combine the following:

    - 10KU lyophilized powder
    - 1000-ul molecular grade TE buffer

    Vortex. Store aliquots @ -20C

4. Lysozyme: CAS 12850-88-3

    To prepare a 10-mg/mL stock, combine the following:
    - 20-mg powder
    - 2000-ul molecular grade TE buffer

    Vortex. Aliquot smaller volumes into tubes. store @ -20C

5. Lysostaphin:

    To prepare an 8-mg/mL (or 4000U/mL) stock, combine the following in the lysostaphin bottle:
    - 2-mg powder (= entire bottle)
    - 250-ul nuclease-free H20

    Vortex. Transfer to a single microfuge tube & store @ -20C

##### (3) Composite Solutions

1. 2X Lysis Solution:

    - 250-ul 0.5M Tris
    - 4-ml 5M NaCl
    - 4-ml 0.5M EDTA
    - 500-ul 20% SLS
    - 1-ml 10% Brij
    - 1-ml 10% Deoxycholate

    The values above are to prepare a total volume of 10-mL. Extrapolate from this formula to prepare your desired volume.


2. 1X Lysis Solution: Perform 1:2 dilution of 2X Lysis solution with H20 [i.e. 25-mL 2X lysis plus 25-mL sterile H20].

    *Note*: Lytic enzymes [Lysozyme, Lysostaphin, etc] plus RNase will be added to both 2X Lysis and 1X Lysis on the day of plug-making (see appendix 1. worksheet)

3. ESP Stock Solution:

    - 10-ml 0.5M TRIS
    - 1-ml 0.5M EDTA

    Add water to 500-mL. Autoclave – store at RT.

    **On day of use**, we need 3-ml ESP solution per plug/sample; so combine the following for each mL total volume:

    - 100-ul 10% SDS
    - 5-ul Proteinase K (20 mg/ml stock)
    - 900-ul ESP stock solution

4. TE Wash Solution:

    - 10-ml 0.5M TRIS
    - 100-ul 0.5M EDTA

    Add water to 500-ml. Autoclave. Store at RT.


5. 1.6% LMP agarose:

    - 1.6g LMP (low melting point) agarose
    - 100-mL TE buffer (10mM tris, 1mM EDTA, pH8.0)

    Microwave for 2-min (stop frequently to swirl solution). Heating step is complete when no agarose particles are seen floating in solution. Temper to 55C before use.

### Sample-type Specific Initial Steps for Tongue Dorsum Samples

* Vortex sample until uniform suspension
* Centrifuge, 10,000 RPM for 15min to pellet microbial cells to bottom of tube.
* Carefully decant/remove supernatant from each of the tubes. Save pellet.

### Procedure

The procedure covers multiple days and distinct steps.

#### Plug Making (Day 1)

{:.notice}
Note: use PFGE plug making worksheet (Appendix 1.) to calculate how to prepare the various solutions based on your desired sample number.

{:.warning}
Need 55C heating block prior to beginning the procedure

1. Melt 1.6% LMP agarose in microwave for 1-min, then temper on the heating block for 3 to 5-min.

2. Pipette 300-ul agarose into 1.5-mL microfuge tubes in 55C heating block.

3. Resuspend a sample (i.e microbial pellet) in 300-ul 2X LYSIS solution (w/ enzymes added).

4. Vortex until a uniform suspension is obtained.

5. Using 1-ml wide-bore pipette tips, add 300-ul cell suspension into 300-ul agarose.

6. Vortex 2-sec to completely mix solutions & return to heat block.

7. Using 1-ml wide-bore pipette tips, use all solution to fill plug mold (~600-ul). Add sample slowly down plug mold side to minimize inclusion of air bubbles.

8. Repeat steps #3 - #7 with any remaining cultures.

9. Refrigerate plug mold for 10-min.

#### Bacterial (Plug) Lysis, part I (Day 1/2)

1. Prepare fresh 1X LYSIS solution with enzymes. Prepare solution according to instructions in PFGE plug making worksheet (Appendix 1.)

2. Add 3-ml 1X LYSIS solution with enzymes to a new set of labeled tubes. If your tubes are labeled with a tape label, transfer of label from the initial culture tube to new culture tube is easy.

3. Remove clips from plug mold.

4. Gently slide mold apart.

5. With flamed spatula, transfer plugs into tubes with 1X LYSIS solution. 14-mL culture tubes work well here.

6. Place tubes on rotator (65-80 RPM), and incubate at 35-37C for minimum 2 hours. This step can proceed overnight, if needed, to accommodate scheduling. Can use lab shaker/incubator for this step.

#### Bacterial (Plug) Lysis, part II (Day 2)

Ensure a 50-55C water bath is available & ready!

1. Prepare ESP solution with enzymes (need 3-mL per sample). Prepare ESP solution according to instructions in PFGE plug making worksheet (Appendix 1.)

3. Decant 1X Lysis solution; leave plug in culture tube.

4. Add 3-ml ESP with enzymes to tube.

5. Place in 50-55C water bath for 1-hr minimum.

6. Remove to room temp for about 5-min before continuing to next step. This time is necessary for the plug to solidify a bit otherwise there is great potential for breakage.

#### TE Wash #1 (Day#2)

1. Decant ESP solution; leave plug in tube.
2. Add ~4-ml TE to culture tube.
3. Place in 50-55C water bath for 1-hr minimum.
4. Remove to room temp for about 5-min before continuing to next step. This time is necessary for the plug to solidify a bit otherwise there is great potential for breakage.

#### F. TE Wash #2 (Day#2)

1. Add ~25-ml of TE to a petri dish (1 dish per sample)

2. Decant TE from culture tube (so that plug remains in tube)

3. Gently knock plug into plastic petri dish.

4. Transfer sample label to dish

5. Place on rotator (65-rpm) at 35C for 30-min minimum.

#### Beta-agarase I digestion

Agarose Digestion and DNA purification will be performed by the [NEB protocol](https://www.neb.com/protocols/0001/01/01/dna-purification-from-agarose-gels-using-beta-agarase-i-m0392) "DNA purification from agarose gels using Beta agarose I (NEB# M0392)".

{:.notice}
Before starting you will need to prepare 100% isopropanol and 70% isopropanol - store both at -20C until needed.

1. With plugs still in the petri dish, aspirate off the TE with a sterile transfer pipette.

2. Using a flamed spatula, cut the plug into 3 fragments (in the petri dish).

3. Transfer the plug pieces to a sterile 1.5-mL microfuge tube.

    - Push them as far down into the tube as possible.
    - When heating, pieces at the very bottom will melt most quickly.

4. Equilibrate the DNA-containing LMP agarose with 2 volumes (~1.0 to 1.5-ml) of 1X Beta-agarase I buffer on ice for 30-min.

    - used 120-ul 10X B-agarase buffer into 1080-ul nuclease-free water.

5. Repeat wash step 1X.

    - used 60-ul 10X B-agarase buffer into 540-ul nuclease-free water

6. Remove the buffer and melt the agarose by incubation at 65C for 10-min.

    - Swirl the tube(s) multiple times to ensure complete melting.
    - Note: Agarases digest denatured/melted agarose, not chunks of agarose

7. Cool fully to 42C.

    - Allow sufficient time for the molten agarose to fully equilibrate as the Beta-agarase enzyme is inactivated at temperatures above 45C.

8. Incubate the molten agarose with 6-ul Beta-agarase I enzyme (i.e. 3 units, 0.5U/ul) at 42C for 1-hour to overnight.

    - Use 1U of agarase per 200-ul gel slice
    - Mix end-over-end rotation 3X

9. Purify and Concentrate the DNA:

    - Add 303-ul 5M ammonium acetate [want 2.5M final concentration] & invert to mix.
    - Using ammonium acetate avoids co-precipitation of oligosaccharides.
    - This adjusts the salt concentration accordingly in preparation for isopropanol precipitation.

    * Chill on ice for 15 min.

    * Centrifuge 10,000 x g for 15-min to pellet undigested carbohydrates, epithelial cells, etc (note: with only salt added, the DNA has not yet precipitated and will be present in the supernatant).

    * Transfer DNA-containing supernatant to new tube with wide-bore pipet tip & very slow pipetting.

        - Split the supernatent into 2 microfuge tubes (500-ul and 250-ul) to accommodate 2 volumes of isopropanol in next steps.
        - Noted white precipitate floating at the liquid surface which did not pellet.

    * Add 2 volumes ice-cold 100% isopropanol – mix gently by end over end rotation.

    * Incubate at -20C for 20-min.

    * Centrifuge 10,000 x g for 15-min to pellet DNA.

    * Remove the supernatant (save pellet).

    * Wash the pellet with 1-mL cold 70% isopropanol. Centrifuge 10,000 x g for 2-min. Remove supernatant.

    * Dry the pellet at room temperature for 5-10 min (keep the tubes in an inverted position). Overdrying the DNA will result in difficulty in resuspending it.

    * Resuspend the sample:
        - If sample was split earlier, resuspend each tube in 50ul elution buffer (EB; 10 mM Tris-HCl pH 8.5) and then recombine into one sample tube. Use wide-bore pipet tips.
        - If sample was NOT split earlier, resuspend the pellet in 100ul elution buffer (EB; 10 mM Tris-HCl pH 8.5)

## Appendix 1. table

{:.notice}
Always prepare extra volumes!

**2X Lysis Solution Component (need 1-mL per sample)**

| Component | Stock Concentration | Final Concentration |Volume required per plug |
| -------- | -------- | -------- |-------- |
| 2X Lysis Solution    | 2X    | 2X    |1-mL    |
| Lysozyme | 10-mg/mL | 1-mg/mL |100-ul |
| RNase    | 10-mg/mL    | 30-mcg/mL     |3-ul    |

**1X Lysis Solution Component (need 3-mL per sample)**

| Component | Stock Concentration | Final Concentration |Volume required per plug |
| -------- | -------- | -------- |-------- |
| 1X Lysis Solution    | 1X    | 1X    |3-mL    |
| Lysozyme | 10-mg/mL | 0.5-mg/mL |150-ul |
| RNase    | 10-mg/mL    | 100-mcg/mL     |10-ul    |
| Lysostaphin    | 8-mg/mL (aka 4000U/mL)   | 50U/mL     |12.5-ul    |
| Mutanolysin    | 10KU/mL    | 0.3KU/mL     |90-ul    |

**ESP Solution Component (3-mL per sample)**

| Component | Stock Concentration | Final Concentration |Volume required per plug |
| -------- | -------- | -------- |-------- |
| ESP sol'n w/o enzymes  | 1X    | 1X    |2.7-mL    |
| 10% SDS | 10% soln | 1% soln |0.3-mL |
|ProteinaseK   | 20-mg/mL    | 100-mcg/mL     |15-ul
