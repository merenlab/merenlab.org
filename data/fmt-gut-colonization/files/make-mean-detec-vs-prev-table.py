<<<<<<< HEAD

#!/usr/bin/env python3
=======
>>>>>>> 05b25422ae07009bc26f00953a8f33577d51e51c
import pandas as pd

donors_meta = pd.read_csv("metadata-donor.txt", sep="\t")

recipients_meta = pd.read_csv("metadata-recipient.txt", sep="\t")

transplants = pd.read_csv("metadata-transplants.txt", sep="\t")

for group in ["DA", "DB"]:

    detection = pd.read_csv(f"detection-FMT-{group}.txt", sep="\t", index_col=0)
    prevalence = pd.read_csv(f"detection-global-by-country-{group}.txt", sep="\t", index_col=0)

    # DONOR

    donor_samples = donors_meta[donors_meta['Group'] == group]
    donor_samples = sorted(set(donor_samples['Sample Name']))
    donor_samples_detection = detection.loc[ : , donor_samples]
    donor_samples_detection["mean_donor_detection"] = donor_samples_detection.mean(axis=1)
    donor_means = donor_samples_detection["mean_donor_detection"]

    # RECIPIENT_POST

    recipient_post_samples = recipients_meta[(recipients_meta['Days Relative to Transplant'] > 0) & \
                                             (recipients_meta['Group'] == group)]
    recipient_post_samples = sorted(set(recipient_post_samples['Sample Name']))
    post_samples_detection = detection.loc[ : , recipient_post_samples]
    post_samples_detection["mean_post_detection"] = post_samples_detection.mean(axis=1)
    post_means = post_samples_detection["mean_post_detection"]

    # MERGING
    df = pd.concat([donor_means, post_means], axis=1)
    final = df.merge(prevalence, how="inner", left_index=True, right_index=True)

    final.to_csv(f"mean-detection-vs-prevalence-{group}.txt", sep = "\t", index=True)
