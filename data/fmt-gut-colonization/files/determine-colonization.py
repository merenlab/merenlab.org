#!/usr/bin/env python3
import pandas as pd
import os.path

detection_cutoff = 0.25     # MAG detection greater than or equal to this value.
coverage_cutoff = 10        # MAG coverage greater than or equal to this value.
abundance_cutoff = 0.01     # MAG abundance greater than or equal to this value.
days_post_cutoff = 7        # POST samples taken greater than this number of days post-FMT.

donors_raw  = pd.read_csv('metadata-donor.txt', sep='\t')
recipients_raw = pd.read_csv('metadata-recipient.txt', sep='\t')
transplants_raw = pd.read_csv('metadata-transplants.txt', sep='\t')

for group in ['DA', 'DB']:
    detection = pd.read_csv(f'detection-FMT-{group}.txt', sep='\t', index_col=0)
    prevalence = pd.read_csv(f'detection-global-by-country-{group}.txt', sep='\t', index_col=0)
    scg_cov = pd.read_csv(f'scg-cov-{group}.txt', sep='\t', index_col=0)
    subpop_num = pd.read_csv(f'subpop-num-{group}.txt', sep='\t', index_col=0)
    desman = pd.read_csv(f'subpop-comp-{group}.txt', sep='\t', na_values='NA')

    donors_meta = donors_raw[donors_raw['Group'] == group]
    recipients_meta = recipients_raw[recipients_raw['Group'] == group]
   
    # Make lists of samples used:
    recipient_samples = sorted(set(recipients_meta['Sample Name']))
    donor_samples = sorted(set(donors_meta['Sample Name']))
    all_samples = donor_samples + recipient_samples
    
    recipients = sorted(set(recipients_meta['Recipient']))
    
    # Make a list of all MAGs:
    MAGs_all = list(detection.index.values.tolist())
    
    # Make a dictionary of PRE samples for each recipient:
    recipient_pre = {}
    for recipient in recipients:
        pre_df = recipients_meta.loc[(recipients_meta['Recipient'] == recipient) & 
                         (recipients_meta['Days Relative to Transplant'] < 0)]
        pre_samples = list(pre_df['Sample Name'])
        recipient_pre[recipient] = pre_samples
   
    # Make a dictionary of POST samples for each recipient:
    recipient_post = {}
    for recipient in recipients:
        post_df = recipients_meta.loc[(recipients_meta['Recipient'] == recipient) & 
                          (recipients_meta['Days Relative to Transplant'] > days_post_cutoff)]
        post_samples = list(post_df['Sample Name'])
        recipient_post[recipient] = post_samples
   
    ####################################
    # THE MESS (FLOWCHART) STARTS HERE #
    ####################################
      
    # TEST 1: the MAG SCGs must have a mean non-outlier cov >= cutoff in at least one donor sample.

    test1_passed = []   # MAGs with coverage above cutoff in any donor sample.

    covered_donor_samples = {}  # For each MAG, the donor samples with sufficient coverage.

    colonized = {}
    did_not_colonize = {}
    ambiguous = {}

    for MAG in MAGs_all:
        # Get coverage in donor samples:
        donor_test = scg_cov.loc[MAG, donor_samples]
        # Make a list of the donor samples with sufficient coverage:
        mag_covered_donor_samples = list(donor_test[donor_test >= coverage_cutoff].index)
        
        if mag_covered_donor_samples != []:
            test1_passed.append(MAG)
            covered_donor_samples[MAG] = mag_covered_donor_samples
            # passed test 1.
        else:
            ambiguous[MAG] = recipients

    for MAG in test1_passed:

        # LIST 1: the subpops present in the donor.
        if subpop_num.loc[MAG, 'Number of subpopulations'] > 1:
            donor_desman = desman[(desman['MAG'] == MAG) & (desman['Sample'].isin(covered_donor_samples[MAG]))]
            subpop_only = donor_desman.iloc[:,2:]
            filtered = subpop_only[subpop_only >= abundance_cutoff]
            filtered.dropna(axis='columns', how='all', inplace=True) # REMOVE THIS LINE?
            list_1 = filtered.columns.tolist()
        else:
            list_1 = ['subpop_1']

        for recipient in recipients:
            # TEST 2: is the MAG detected in the transplant sample?
            transplant_sample = transplants_raw.loc[(transplants_raw['Group'] == group) &
                                     (transplants_raw['Recipient'] == recipient), 'Sample Name']                                                                                                        
            transplant_sample = transplant_sample.to_string(index=False).strip()
            transplant_detec = detection.loc[MAG, transplant_sample]
            if transplant_detec >= detection_cutoff:
                # passed test 2.                

                # Is the MAG detected post-FMT?
                post_detec = detection.loc[MAG, recipient_post[recipient]]
                if (post_detec >= detection_cutoff).any():     
                    
                    # TEST 3: Is the MAG sufficiently covered in any post-FMT sample?
                    post_cov = scg_cov.loc[MAG, recipient_post[recipient]] 
                    if (post_cov >= coverage_cutoff).any():
                        # passed test 3.

                        # LIST 2: the subpops present in the recipient post-FMT.
                        covered_post_samples = post_cov[post_cov >= coverage_cutoff].index.tolist()
                        if subpop_num.loc[MAG, 'Number of subpopulations'] > 1:
                            post_desman = desman[(desman['MAG'] == MAG) & (desman['Sample'].isin(covered_post_samples))]
                            subpop_only = post_desman.iloc[:,2:]
                            filtered = subpop_only[subpop_only >= abundance_cutoff]
                            filtered.dropna(axis='columns', how='all', inplace=True) # REMOVE THIS LINE?
                            list_2 = filtered.columns.tolist()
                        else:
                            list_2 = ['subpop_1']
                        
                        # Are any donor subpops in the recipient post-FMT?
                        if bool(set(list_1) & set(list_2)):
                            
                            # LIST 3: subpops in donor and recipient post-FMT.
                            list_3 = (set(list_1) & set(list_2))
            
                            # Is the MAG detected pre-FMT?
                            pre_detec = detection.loc[MAG, recipient_pre[recipient]]
                            if (pre_detec >= detection_cutoff).any():
                                
                                # TEST 4: Is the MAG sufficiently covered in any pre-FMT sample?
                                pre_cov = scg_cov.loc[MAG, recipient_pre[recipient]]
                                if (pre_cov >= coverage_cutoff).any():
                                    # passed test 4.

                                   # LIST 4: the subpops present in the recipient pre-FMT.
                                    covered_pre_samples = pre_cov[pre_cov >= coverage_cutoff].index.tolist()
                                    if subpop_num.loc[MAG, 'Number of subpopulations'] > 1:
                                        pre_desman = desman[(desman['MAG'] == MAG) & (desman['Sample'].isin(covered_pre_samples))]
                                        subpop_only = pre_desman.iloc[:,2:]
                                        filtered = subpop_only[subpop_only >= abundance_cutoff]
                                        filtered.dropna(axis='columns', how='all', inplace=True) # REMOVE THIS LINE?
                                        list_4 = filtered.columns.tolist()
                                    else:
                                        list_4 = ['subpop_1']                                    
                                    
                                    # TEST 5: are there any donor subpops in the recipient post-FMT
                                    # that aren't also in the recipient pre-FMT?
                                    if bool(set(list_3) - set(list_4)):
                                        colonized.setdefault(MAG, []).append(recipient)
                                        # passed test 5 (and colonized!)

                                        # NOW REWIND #

                                    else:
                                        # failed test 5.
                                        ambiguous.setdefault(MAG, []).append(recipient)

                                else:
                                    # failed test 4.
                                    ambiguous.setdefault(MAG, []).append(recipient)

                            else:
                                # colonized!
                                colonized.setdefault(MAG, []).append(recipient)

                        else:
                            # did not colonize!
                            did_not_colonize.setdefault(MAG, []).append(recipient)

                    else:
                        # failed test 3.
                        ambiguous.setdefault(MAG, []).append(recipient)
                    
                else:
                    # did not colonize!
                    did_not_colonize.setdefault(MAG, []).append(recipient)
            
            else:
                # failed test 2.
                ambiguous.setdefault(MAG, []).append(recipient)


    # Write simple results to file:
    file_colo = f'colonized-{group}.txt' 
    file_no_colo = f'did-not-colonize-{group}.txt'
    
    with open(file_colo, 'w') as filetowrite:
        filetowrite.write('MAG\trecipient\n')
        for MAG in colonized:
            for recipient in colonized[MAG]:
                filetowrite.write(MAG + '\t' + recipient + '\n')

    with open(file_no_colo, 'w') as filetowrite:
        filetowrite.write('MAG\trecipient\n')
        for MAG in did_not_colonize:
            for recipient in did_not_colonize[MAG]:
                filetowrite.write(MAG + '\t' + recipient + '\n')
