#!/usr/bin/python
"""Functions to help with execution of system commands."""

import os
import sys
import json
import shutil
import logging
import traceback
import subprocess


def run_cmds(commands, retry=0, catchExcept=False, stdout=None):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    if stdout is None:
        p = subprocess.Popen(commands,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        stdout, stderr = p.communicate()
    else:
        with open(stdout, "wt") as fo:
            p = subprocess.Popen(commands,
                                 stderr=subprocess.PIPE,
                                 stdout=fo)
            stdout, stderr = p.communicate()
        stdout = False
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.split('\n'):
            logging.info(line)

    # Check the exit code
    if exitcode != 0 and retry > 0:
        msg = "Exit code {}, retrying {} more times".format(exitcode, retry)
        logging.info(msg)
        run_cmds(commands, retry=retry - 1)
    elif exitcode != 0 and catchExcept:
        msg = "Exit code was {}, but we will continue anyway"
        logging.info(msg.format(exitcode))
    else:
        assert exitcode == 0, "Exit code {}".format(exitcode)


def align_reads(read_fp,               # FASTQ file path
                db_fp,                 # Local path to DB
                temp_folder,           # Folder for results
                query_gencode=11,      # Genetic code
                threads=1,             # Threads
                blocks=4):             # Memory block size

    """Align a set of reads with DIAMOND."""

    align_fp = "{}.sam".format(read_fp)
    logging.info("Input reads: {}".format(read_fp))
    logging.info("Reference database: {}".format(db_fp))
    logging.info("Genetic code: {}".format(query_gencode))
    logging.info("Threads: {}".format(threads))
    logging.info("Output: {}".format(align_fp))

    run_cmds([
            "diamond",
            "blastx",
            "--query", read_fp,             # Input FASTQ
            "--out", align_fp,              # Alignment file
            "--threads", str(threads),      # Threads
            "--db", db_fp,                  # Reference database
            "--outfmt", "6",                # Output format
            "qseqid", "sseqid",
            "pident", "length",
            "qstart", "qend",
            "sstart", "send",
            "evalue", "bitscore",
            "qlen", "slen",
            "--top", "10",                  # Alignments within 10% of the best
            "--block-size", str(blocks),    # Memory block size
            "--query-gencode",              # Genetic code
            str(query_gencode),
            "--unal", "0",                  # Don't report unaligned reads
            ])

    return align_fp


def get_reference_database(ref_db, temp_folder, ending=None):
    """Get a reference database folder."""
    assert ref_db is not None, "Must provide reference database path"
    assert temp_folder is not None, "Must provide temp folder"
    if ending is not None:
        msg = "Ref DB must end with {} ({})"
        assert ref_db.endswith(ending), msg.format(ending, ref_db)

    # Get files from AWS S3
    if ref_db.startswith('s3://'):
        logging.info("Getting reference database from S3: " + ref_db)

        # Save the database to the local temp folder
        local_fp = os.path.join(
            temp_folder,
            ref_db.split('/')[-1]
        )

        assert os.path.exists(local_fp) is False

        logging.info("Saving database to " + local_fp)
        run_cmds([
            'aws',
            's3',
            'cp',
            '--quiet',
            '--sse',
            'AES256',
            ref_db,
            local_fp
        ])

        return local_fp

    else:
        # Treat the input as a local path
        logging.info("Getting reference database from local path: " + ref_db)

        assert os.path.exists(ref_db)

        return ref_db


def return_alignments(align_fp, output_path):
    """Return the alignment file."""
    logging.info("Compressing " + align_fp)
    run_cmds([
        "pigz", align_fp
    ])
    align_fp = align_fp + ".gz"
    assert os.path.exists(align_fp)
    
    # Make sure that the file has a consistent ending
    assert align_fp.endswith(".sam.gz")
    assert output_path.endswith(".sam.gz")
    
    if output_path.startswith("s3:/"):
        # Copy to S3
        run_cmds([
            'aws',
            's3',
            'cp',
            '--quiet',
            '--sse',
            'AES256',
            align_fp,
            output_path
        ])
        os.unlink(align_fp)
    else:
        # Copy to local folder
        run_cmds(['mv', align_fp, output_path])


def return_results(out, output_path, temp_folder):
    """Write out the final results as a JSON object."""
    # Make a temporary file
    temp_fp = os.path.join(temp_folder, 'temp.json')
    with open(temp_fp, 'wt') as fo:
        json.dump(out, fo)

    # Compress the output
    run_cmds(['gzip', temp_fp])
    temp_fp = temp_fp + '.gz'

    if output_path.startswith('s3://'):
        # Copy to S3
        run_cmds([
            'aws',
            's3',
            'cp',
            '--quiet',
            '--sse',
            'AES256',
            temp_fp,
            output_path])
        os.unlink(temp_fp)
    else:
        # Copy to local folder
        run_cmds(['mv', temp_fp, output_path])


def exit_and_clean_up(temp_folder):
    """Log the error messages and delete the temporary folder."""
    # Capture the traceback
    logging.info("There was an unexpected failure")
    exc_type, exc_value, exc_traceback = sys.exc_info()
    for line in traceback.format_tb(exc_traceback):
        logging.info(line)

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Exit
    logging.info("Exit type: {}".format(exc_type))
    logging.info("Exit code: {}".format(exc_value))
    sys.exit(exc_value)
