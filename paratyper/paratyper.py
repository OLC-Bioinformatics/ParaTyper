#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, SetupLogging, \
    run_subprocess, write_to_logfile
from genemethods.sipprCommon.objectprep import Objectprep
from argparse import ArgumentParser
from threading import Thread
from queue import Queue
import multiprocessing
import json
import time
import os
__author__ = 'adamkoziol'


class ParasiteTyping(object):

    def main(self):
        self.object_prep()
        self.normalise_reads()
        self.sketch_reads()
        self.mash_reads()

    def object_prep(self):
        """
        Create metadata objects for every .ab1 file in the supplied sequence path
        """
        # Create the objects to be used in the analyses
        objects = Objectprep(self)
        objects.objectprep()
        self.metadata = [sample for sample in objects.samples.samples]

    def normalise_reads(self):
        """

        """
        for sample in self.metadata:
            sample.general.normalised_reads = os.path.join(sample.general.outputdirectory, '{sn}_normalised.fastq.gz'
                                                           .format(sn=sample.name))
            sample.commands.normalise = 'bbnorm.sh in1={forward} in2={reverse} maxdepth=50 mindepth=10 ecc=t out={out}'\
                .format(forward=sample.general.fastqfiles[0],
                        reverse=sample.general.fastqfiles[1],
                        out=sample.general.normalised_reads)
            if not os.path.isfile(sample.general.normalised_reads):
                out, err = run_subprocess(sample.commands.normalise)
                write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample.commands.normalise,
                                                           out=out),
                                 err=err,
                                 logfile=self.logfile,
                                 samplelog=sample.general.logout,
                                 sampleerr=sample.general.logerr,
                                 analysislog=None,
                                 analysiserr=None)

    def sketch_reads(self):
        """

        """
        # Create the threads for the analysis
        for i in range(self.cpus):
            threads = Thread(target=self.sketch, args=())
            threads.setDaemon(True)
            threads.start()
        for sample in self.metadata:
            # Create the analysis type-specific GenObject
            setattr(sample, 'paratyper', GenObject())
            sample.paratyper.sketchfilenoext = os.path.join(sample.general.outputdirectory, sample.name)
            sample.paratyper.sketchfile = sample.paratyper.sketchfilenoext + '.msh'
            sample.commands.sketch = 'mash sketch -m 2 -r {reads} -o {output_file}' \
                .format(reads=sample.general.normalised_reads,
                        output_file=sample.paratyper.sketchfilenoext)
            self.sketchqueue.put(sample)
        # Join the threads
        self.sketchqueue.join()

    def sketch(self):
        while True:
            sample = self.sketchqueue.get()
            if not os.path.isfile(sample.paratyper.sketchfile):
                # Run the command
                out, err = run_subprocess(sample.commands.sketch)
                write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample.commands.sketch,
                                                           out=out),
                                 err=err,
                                 logfile=self.logfile,
                                 samplelog=sample.general.logout,
                                 sampleerr=sample.general.logerr,
                                 analysislog=None,
                                 analysiserr=None)
            self.sketchqueue.task_done()

    def mash_reads(self):
        # Create the threads for the analysis
        for i in range(self.cpus):
                threads = Thread(target=self.mash, args=())
                threads.setDaemon(True)
                threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            sample.paratyper.mashresults = os.path.join(sample.general.outputdirectory, '{sn}.tab'
                                                        .format(sn=sample.name))
            sample.commands.mash = \
                'mash dist {ref_sketch} {sample_sketch} | sort -gk3 > {results}'\
                .format(ref_sketch=self.reference_mash_sketch_file,
                        sample_sketch=sample.paratyper.sketchfile,
                        results=sample.paratyper.mashresults)
            self.mashqueue.put(sample)
        # Join the threads
        self.mashqueue.join()

    def mash(self):
        while True:
            sample = self.mashqueue.get()
            #
            if not os.path.isfile(sample.paratyper.mashresults):
                # Run the command
                out, err = run_subprocess(sample.commands.mash)
                write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample.commands.mash,
                                                           out=out),
                                 err=err,
                                 logfile=self.logfile,
                                 samplelog=sample.general.logout,
                                 sampleerr=sample.general.logerr,
                                 analysislog=None,
                                 analysiserr=None)
            self.mashqueue.task_done()

    def parse_mash_outputs(self):
        """

        """
        pass

    def report(self):
        """
        Create a JSON-formatted report
        """
        for sample in self.metadata:
            # Add the sample to the output dictionary as sample name: attribute name: attribute: value
            self.output_dict[sample.name] = sample.dump()
            # Remove the 'unwanted keys' key from the dictionary, as this is only useful for metadata objects
            self.output_dict[sample.name].pop('unwanted_keys', None)
        # Open the metadata file to write
        with open(self.json_report, 'w') as metadatafile:
            # Write the json dump of the object dump to the metadata file
            json.dump(self.output_dict, metadatafile, sort_keys=True, indent=4, separators=(',', ': '))

    def __init__(self, sequencepath, reportpath):
        # Allow for relative paths
        if sequencepath.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(sequencepath)))
        else:
            self.path = os.path.abspath(os.path.join(sequencepath))
        assert os.path.isdir(self.path), 'Cannot locate supplied sequence path: {seq_path}' \
            .format(seq_path=self.path)
        self.sequencepath = self.path
        if reportpath.startswith('~'):
            self.reportpath = os.path.abspath(os.path.expanduser(os.path.join(reportpath)))
        else:
            self.reportpath = os.path.abspath(os.path.join(reportpath))
        make_path(self.reportpath)
        assert os.path.isdir(self.reportpath), 'Could not create the requested report directory: {rep_path}' \
            .format(rep_path=self.reportpath)
        # Define the start time for legacy code compatibility
        self.starttime = time.time()
        self.logfile = os.path.join(self.path, 'log')
        self.cpus = multiprocessing.cpu_count() - 1
        self.sketchqueue = Queue(maxsize=self.cpus)
        self.mashqueue = Queue(maxsize=self.cpus)
        # Extract the path of this file - will be used to find the necessary accessory files
        self.homepath = os.path.split(os.path.abspath(__file__))[0]
        # self.reference_mash_sketch_file = os.path.join(self.homepath, 'toxoplasma.msh')
        self.reference_mash_sketch_file = '/mnt/nas2/redmine/bio_requests/17874/reference_sequences/toxoplasma.msh'
        self.metadata = list()
        self.output_dict = dict()
        self.json_report = os.path.join(self.reportpath, 'para_typer_outputs.json')


def main():
    parser = ArgumentParser(description='Perform parasite typing')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path of folder containing FASTQ files to process.')
    parser.add_argument('-r', '--reportpath',
                        required=True,
                        help='Path in which reports are to be created')

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Allow debug-level logging to be printed to the terminal')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=arguments.verbose)
    para_typer = ParasiteTyping(sequencepath=arguments.sequencepath,
                                reportpath=arguments.reportpath)
    para_typer.main()


if __name__ == '__main__':
    main()
