#!/usr/bin/env bats

@test "SRA Toolkit v2.8.2" {
  v="$(fastq-dump --version)"
  [[ "$v" =~ "2.8.2" ]]
}


@test "AWS CLI v1.11.146" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "1.11.146" ]]
}


@test "Curl v7.47.0" {
  v="$(curl --version)"
  [[ "$v" =~ "7.47.0" ]]
}

@test "DIAMOND v0.9.10" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.10" ]]
}

@test "Make sure the run script is in the PATH" {
  h="$(map_viruses.py -h 2>&1)"

  [[ "$h" =~ "Align a set of reads" ]]
}

@test "Alignment parser" {
  h="$(python /usr/map_viruses/lib/test_alignment_parser.py)"

  [[ "$h" =~ "Success" ]]
}

