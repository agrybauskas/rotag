#!/bin/bash
cd "$(dirname "$0")"

angle_range="0-2"
angle_df=0.2

../programs/sample_angles "${angle_range}" "${angle_df}" | sort -n
