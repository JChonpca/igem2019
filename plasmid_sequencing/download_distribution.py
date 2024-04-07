import os
import csv
import time
import subprocess

input_folder="plates"
output_folder="output"

# Iterate through files in the input folder
for filename in os.listdir(input_folder):
    filepath = os.path.join(input_folder, filename)
    
    # Check if the path points to a file and ends with '.csv'
    if os.path.isfile(filepath) and filename.lower().endswith(".csv"):
        print(f"Reading CSV file: {filename}")

        # Open the CSV file and perform actions (example: print each row)
        with open(filepath, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            header = next(csv_reader)
            for row in csv_reader:
                print(row)
                row[2] = row[2].strip()
                output_file_path = os.path.join(output_folder, row[2]) + ".txt"
                if os.path.exists(output_file_path):
                    continue
                command = "curl -o {} \"https://parts.igem.org/cgi/partsdb/composite_edit/putseq.cgi?part={}\"".format(output_file_path, row[2])
                print("Downloading {} {} {}".format(row[0], row[1], row[2]))
                subprocess.run(command, shell=True, check=True)
                time.sleep(5)

f = open("2018_Distribution.fasta", "w")

# Create merged FASTA with useful names
for filename in os.listdir(input_folder):
    filepath = os.path.join(input_folder, filename)
    
    # Check if the path points to a file and ends with '.csv'
    if os.path.isfile(filepath) and filename.lower().endswith(".csv"):
        print(f"Reading CSV file: {filename}")
        plate_number = filename[15:16]
        print("Plate: " + plate_number)

        # Open the CSV file and perform actions (example: print each row)
        with open(filepath, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            header = next(csv_reader)
            for row in csv_reader:
                print(row)

                row[2] = row[2].strip()
                row[1] = row[1].strip()

                output_file_path = os.path.join(output_folder, row[2]) + ".txt"
                lines = []
                with open(output_file_path) as file:
                    lines = [line.rstrip() for line in file]

                f.write(">{}_P{}_{}\n".format(row[2], plate_number, row[1]))
                f.write(lines[1] + "\n") 



