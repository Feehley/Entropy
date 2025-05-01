#!/usr/bin/env python3

from math import log

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from argparse import ArgumentParser
from io import BytesIO as IO
from os import get_terminal_size, path, getcwd
from numpy import ceil, frombuffer, unpackbits

from PIL import Image


class Shannon_Entropy():
    def __init__(self, **kwargs):
        if len(kwargs) > 1:
            self.filename = kwargs["input_file"]
            self.xxd      = kwargs["xxd"]
            self.graph    = kwargs["graph"]
            self.raster   = kwargs["raster"]
            self.compfile = kwargs["compfile"]
            self.output   = kwargs["output"]
            self.clio     = kwargs["clio"]
        else:
            self.xxd      = False
            self.graph    = False
            self.raster   = False
            self.compfile = None
            self.output   = False
            self.clio     = False

        self.data = None
        self.data_bytes = None
        self.entropy_score = 0.0
        self.normalized = []
        self.ylim_norm = ((32 - 0)/(256 - 0))
        self.ymax_norm = ((127 - 0)/(256 - 0))
        self.cutline = "#" * get_terminal_size()[0]
        self.frame_width = 1000

        if len(kwargs) > 1:
            self.check_paths()
            self.check_xxd()
            self.check_graph()
            self.check_raster()
            self.check_compfile()


    def check_paths(self):
        if path.exists(self.filename):
            self.data = self.read_file(self.filename)
            self.entropy_score = self.shannon()
            print(self.cutline)
            print(f"[+] {self.filename}:\n\tShannon Entropy: {self.entropy_score}")
            print(self.cutline)



    def check_xxd(self):
        if self.xxd:
            self.xxd = self.bytes_to_c_arr()
            print(f"{self.cutline}\n{self.xxd}\n{self.cutline}")


    def check_output(self):
        if self.output:
            filename = self.clean_filename(self.filename)
            filename = f"{filename}.xxd"
            with open(filename, 'wb') as fptr:
                fptr.write(f"{self.xxd}\n")


    def check_graph(self):
        if self.graph:
            self.normalize()
            self.plot_entropy()


    def check_raster(self):
        if self.raster:
            self.generate_bit_raster()


    def check_compfile(self):
        if self.compfile != None:
            if path.exists(self.compfile):
                self.comp_data = self.read_file(self.compfile)
                self.compare_entropy()
            else:
                print("[!] Please enter a real file path.")


    # Clean file names
    @staticmethod
    # self.clean_filename(
    def clean_filename(filename):
        if "/" in filename: filename = filename.split("/")[-1]
        if "." in filename: filename = filename.split(".")[0]
        return filename


    # Read in a file and convert it to an array
    def read_file(self, filename):
        with open(filename, 'rb') as fptr:
            data = fptr.read()
        return data


    # Converting every piece of data into something between 0 and 1
    def normalize(self):
        xmin = 0
        xmax = 256

        for dat in self.data:
            if type(dat) == str:
                self.normalized.append((ord(dat) - xmin)/(xmax - xmin))
            else:
                self.normalized.append((dat - xmin)/(xmax - xmin))

        return self.normalized


    # Creating hex repr of ticker data
    @staticmethod
    def to_hex(x, pop):
        return "%x" % self.data


    # Creating byte values
    def to_byte_array(self, data):
        if type(data) != bytes:
            io = IO(bytes(data, "utf-8")).read(len(data))
            fb = unpackbits(frombuffer(io, dtype="uint8"))
        else:
            io = IO(data).read(len(data))
            fb = unpackbits(frombuffer(io, dtype="uint8"))

        data_bytes = fb
        return data_bytes


    # Performs Shannon Entropy analysis on a given block of data
    def shannon(self):
        if self.data:
            length = len(self.data)

            seen = dict(((chr(x), 0) for x in range(0, 256)))
            for byte in self.data:
                if type(byte) != str:
                    byte = chr(byte)
                seen[byte] += 1

            for x in range(0, 256):
                p_x = float(seen[chr(x)]) / length
                if p_x > 0: self.entropy_score -= p_x * log(p_x, 2)

        self.entropy_score = self.entropy_score / 8
        return self.entropy_score


    # Creating graph
    def plot_entropy(self):
        fmt = ticker.FuncFormatter(self.to_hex)

        plt.figure(figsize=(12, 6))
        plt.ylim([0, 1])

        ax = plt.gca()
        ax.set_title(self.filename)
        plt.axhspan(self.ylim_norm, self.ymax_norm, facecolor="green", alpha=0.5)
        plt.plot(self.normalized)

        xlabels = map(lambda t: "0x%08X" % int(t), ax.get_xticks())
        ax.set_xticks([x for x in range(0, len(list(xlabels)))])
        ax.set_xticklabels(xlabels)

        if self.output:
            filename = self.clean_filename(self.filename)
            filename = f"{getcwd()}/{filename}_entropy.png"
            plt.savefig(filename)
        plt.show()
        return


    # Bitraster
    def generate_bit_raster(self):
        self.data_bytes = self.to_byte_array(self.data)
        num_frames = int(ceil(len(self.data_bytes)/self.frame_width))
        zero = 0
        one  = 1

        img = Image.new("1", (self.frame_width, num_frames), color = zero)
        img.putdata(self.data_bytes, 1, 0)

        if self.clio:
            term_width = 80
            zero = u'\u2591'
            one  = u'\u2593'
            tmp  = []

            for i in range(0, 3999, term_width):
                line = self.data_bytes[i:i+term_width]
                l_no = f"{str(hex(i+term_width)).split('x')[1].zfill(4)}"
                tmp  = [l_no, "\t"]

                for j in range(len(line)):
                    if line[j] == 0:
                        tmp.append(zero)
                    elif line[j] == 1:
                        tmp.append(one)

                tmp = [str(t) for t in tmp]
                print(''.join(tmp))
        else:
            img.show()

        if self.output:
            filename = self.clean_filename(self.filename)
            filename = f"{getcwd()}/{filename}_bitraster.png"
            img.save(filename, "PNG")
        return


    # Compare Bitrasters
    def compare_entropy(self):
        base_file_bytes    = self.to_byte_array(self.data)
        compare_file_bytes = self.to_byte_array(self.comp_data)
        tmp_file_bytes     = []

        for i, j in zip(base_file_bytes, compare_file_bytes):
            if i != j:
                if i > j:
                    tmp_file_bytes.append(bytes(i - j))
                else:
                    tmp_file_bytes.append(bytes(j - i))
            else:
                tmp_file_bytes.append(b'0')

        tmp_file_bytes = b''.join(tmp_file_bytes)
        num_frames = int(ceil(len(tmp_file_bytes)/self.frame_width))
        zero = 0
        one  = 1

        img = Image.new('1', (self.frame_width, num_frames), color = zero)
        img.putdata(tmp_file_bytes, 1, 0)

        if self.clio:
            term_width = 80
            zero = u'\u2591'
            one  = u'\u2593'
            tmp  = []

            for i in range(0, 1999, term_width):
                line = tmp_file_bytes[i:i+term_width]
                l_no = f"{str(hex(i+term_width)).split('x')[1].zfill(4)}"
                tmp  = [l_no, "\t"]

                for j in range(len(line)):
                    if line[j] == 0:
                        tmp.append(zero)
                    elif line[j] == 1:
                        tmp.append(one)

                tmp = [str(t) for t in tmp]
                print(''.join(tmp))
        else:
            img.show()

        if self.output:
            filename = self.clean_filename(self.filename)
            compfile = self.clean_filename(self.comp_data)
            filename = f"{getcwd()}/{filename}_vs_{compfile}_bitraster.png"
            img.save(filename, "PNG")


    # XXD look-a-like
    def bytes_to_c_arr(self, lowercase = True):
        hex_arr = []
        data = [format(b, "#04x" if lowercase else "#04X") for b in self.data]

        for i in range(0, len(data), 16):
            start_hex   = str(hex(i)[2:]).zfill(8)
            content_hex = ' '.join(data[i:i+16]).replace("0x", " ")
            tmp         = ""

            for x in self.data[i:i+16]:
                # Printable ASCII Range
                if x <= 127 and x >= 32:
                    tmp += chr(x)
                else:
                    tmp += "."
            hex_arr.append(f"{start_hex}: {content_hex}\t{tmp}")
        return "\n".join(hex_arr)


# Get CLI Input
def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file",
                        help = "Input file")
    parser.add_argument("-o",
                        "--output",
                        action  = "store_true",
                        help    = "Save output(s)",
                        default = False)
    parser.add_argument("-l",
                        "--clio",
                        action  = "store_true",
                        help    = "Preview first 2000 bytes",
                        default = False)
    parser.add_argument("-g",
                        "--graph",
                        action  = "store_true",
                        help    = "Show graphical represenation",
                        default = False)
    parser.add_argument("-r",
                        "--raster",
                        action  = "store_true",
                        help    = "Show bitraster representation",
                        default = False)
    parser.add_argument("-c",
                        "--compfile",
                        action   = "store",
                        help     = "Comparison file; outputs a bitraster",
                        required = False)
    parser.add_argument("-x",
                        "--xxd",
                        action  = "store_true",
                        help    = "Show hex representation",
                        default = False)
    return parser.parse_args()


def main():
    args = get_args()
    Shannon_Entropy(input_file = args.input_file,
                    output     = args.output,
                    clio       = args.clio,
                    graph      = args.graph,
                    raster     = args.raster,
                    compfile   = args.compfile,
                    xxd        = args.xxd)


if __name__ == '__main__':
    exit(main())
