#!/bin/bash
cd "$(dirname "$0")"

angle_range="0-2"
small_angle=0.2

../programs/sample_angles "${angle_range}" "${small_angle}"
