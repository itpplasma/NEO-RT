#!/bin/bash

# shellcheck disable=SC2029

set -euo pipefail

run_tag="${1:?usage: stage_gpec_kinetic.sh RUN_TAG}"
remote="scluster"
remote_root="/home/ert/codex-runs/${run_tag}"
source_root="${remote_root}/source"
repository_root="$(git rev-parse --show-toplevel)"
ideal_root="/mnt/storage/codex-gpec/tc24-gpec-jxb-a430ad2-20260730/baseline"
profile_root="/mnt/storage/codex-gpec/tc24-pentrc-supplier-profilefix-l0-l5-20260726-2792ebd0-2e031ae6-r1/common"
tc24_root="/home/ert/proj/iter_tc24"
gpec_commit="2e031ae675205ec9d69bb0298855155737f72579"

test -f "${profile_root}/profiles.kin"

ssh "${remote}" "test ! -e '${remote_root}' && mkdir -p \
'${remote_root}/debs' '${remote_root}/deps' \
'${remote_root}/input-ideal' '${remote_root}/input-kinetic' \
'${remote_root}/run-ideal' '${remote_root}/run-kinetic'"
ssh "${remote}" "git clone --no-checkout \
https://github.com/PrincetonUniversity/GPEC.git '${source_root}'"
ssh "${remote}" "git -C '${source_root}' fetch origin \
refs/pull/281/head && git -C '${source_root}' checkout --detach '${gpec_commit}'"
ssh "${remote}" \
    "test \"\$(git -C '${source_root}' rev-parse HEAD)\" = '${gpec_commit}'"

for case_name in ideal kinetic; do
    input_root="${remote_root}/input-${case_name}"
    rsync -a \
        "${ideal_root}/coil.in" \
        "${ideal_root}/dcon.in" \
        "${ideal_root}/equil.in" \
        "${ideal_root}/gpec.in" \
        "${ideal_root}/pentrc.in" \
        "${ideal_root}/vac.in" \
        "${remote}:${input_root}/"
    rsync -a \
        "${profile_root}/profiles.kin" \
        "${profile_root}/profile_provenance.json" \
        "${remote}:${input_root}/"
    rsync -a \
        "${tc24_root}/gfiles/ngfile_p1_run_eq_modx03_20260527" \
        "${remote}:${input_root}/"
    rsync -a \
        "${tc24_root}/COILS/" \
        "${remote}:${input_root}/COILS/"

    ssh "${remote}" "perl -0pi -e \
    's/parallel_threads = 32/parallel_threads = 16/' \
    '${input_root}/dcon.in'"
    ssh "${remote}" "perl -0pi -e \
    \"s#eq_filename = .*#eq_filename = \
'${input_root}/ngfile_p1_run_eq_modx03_20260527'#\" \
    '${input_root}/equil.in'"
    ssh "${remote}" "perl -0pi -e \
    \"s#data_dir = .*#data_dir = '${input_root}/COILS'#\" \
    '${input_root}/coil.in'"
done

ssh "${remote}" "env --chdir='${remote_root}/debs' apt-get download \
libnetcdff-dev libnetcdff7"
ssh "${remote}" "for package in '${remote_root}'/debs/*.deb; do \
dpkg-deb -x \"\$package\" '${remote_root}/deps'; done"

ssh "${remote}" "perl -0pi -e \
's/kin_flag = \\.false\\./kin_flag = .true./; \
s/con_flag = \\.false\\./con_flag = .true./' \
'${remote_root}/input-kinetic/dcon.in'"

rsync -a \
    "${repository_root}/MARS-JXB/gpec_kinetic_worker.slurm" \
    "${remote}:${remote_root}/"

ssh "${remote}" "sbatch --parsable \
--output='${remote_root}/slurm-%j.out' \
--export=ALL,CAMPAIGN_ROOT='${remote_root}' \
'${remote_root}/gpec_kinetic_worker.slurm'"
