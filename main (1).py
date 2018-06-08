from pyspark import SparkContext, SparkConf, StorageLevel
from operator import add
import swiftclient.client
import sys
import pysam
import urllib
import config
import os
import time
import hashlib


def downloadFile(filename, remote_path, local_path):
    for i in range(0,10):
    try:
        print "Downloading ", filename
        urllib.urlretrieve(remote_path,local_path)
        downloaded_check_sum = hashlib.md5(open(local_path, 'rb').read()).hexdigest()
        
        if str(check_sum) != str(downloaded_check_sum):
            print "Checksum failed on try ", i, filename
    else:
        print "Checksum succeeded ", filename
        break
    except:
       print "Failed to Download on try ", i, filename


def removeFile(filename):
    extensions = ["", "?", "?.csi"]
    for extension in extensions:
        try:
            os.remove(local_path + extension)
        except OSError:
            pass


def process(bam_file):
    filename  = bam_file['name']
    check_sum = bam_file['hash']
        
    remote_path = "http://130.238.29.253:8080/swift/v1/1000-genomes-dataset/"+filename
    local_path = "/home/ubuntu/" + filename  
    
    start_time_download = time.time()  
    downloadFile(filename, remote_path, local_path)

    start_time_mapping = time.time()

    result_list = []
    K = 10

    # Process the downloaded file
    with pysam.AlignmentFile(local_path,"rb") as samfile:
        j = 0
        try:
            data = samfile.fetch(until_eof=True)
            for r in data:
                j+=1

                # Due to a memory error, the number of lines in a single
                # file is capped to 1000000
                if j > 1000000:
                    break

                if r.is_unmapped and not(r.mate_is_unmapped):
                    result_list.append(('POSITION', (r.reference_start, 1)))
                    start_pos = r.query_alignment_start
                    end_pos = r.query_alignment_end

                    # Produce k-mers
                    for i in range(start_pos, end_pos - (K-1)):
                        sequence = r.query_sequence[i:i+K]
                        if 'N' not in sequence:
                            result_list.append(('KMER', (sequence, 1)))
        except IOError as e:
            print e
            pass

    end_time = time.time()
    result_list.append(('TIME-DOWNLOAD', (1, end_time - start_time_download)))
    result_list.append(('TIME-MAPPING', (1, end_time - start_time_mapping)))

    removeFile(filename)

    return result_list


def main():
    # Initialize Spark
    configuration = SparkConf().setAppName("1000-genomes Project")
    spark_context = SparkContext(conf=configuration)

    container_name = "1000-genomes-dataset"
    config = {'user':os.environ['OS_USERNAME'], 
          'key':os.environ['OS_PASSWORD'],
          'tenant_name':os.environ['OS_TENANT_NAME'],
          'authurl':os.environ['OS_AUTH_URL']}

    # Connect to Object-Storage to retrieve filenames
    conn = swiftclient.client.Connection(auth_version=3, **config)
    (storage_url, auth_token) = conn.get_auth()
    (response, content) = swiftclient.client.get_container(url=storage_url,container=container_name, token=auth_token)
    
    num_files_to_include = 3
    names = filter(lambda t: t['name'][-4:] == '.bam', content)   
    filelist = [{"name" : c['name'].strip(), "hash" : c['hash'].strip()} for c in names[:num_files_to_include]]
    numpart = len(filelist)/2

    # Distribute the processing of the files
    filenames = spark_context.parallelize(filelist, numpart)
    mapped_data = filenames.flatMap(process).persist(StorageLevel(True, True, False, True, 1))

    # Filter the different results returned from the processing 
    start_time_filtering = time.time()
    kmers = mapped_data.filter(lambda (k, (v, e)): k == "KMER").map(lambda (k, v): v).reduceByKey(add,numPartitions=5)
    positions = mapped_data.filter(lambda (k,v): k == "POSITION").map(lambda (k, v): v).reduceByKey(add,numPartitions=5) 

    time_filtering = time.time() - start_time_filtering
    time_mapping = mapped_data.filter(lambda (k,v): k == "TIME-MAPPING").map(lambda (k, v): v).reduceByKey(add,numPartitions=5)
    time_download = mapped_data.filter(lambda (k,v): k == "TIME-DOWNLOAD").map(lambda (k, v): v).reduceByKey(add,numPartitions=5)

    time_mapping = time_mapping.collect()
    time_download = time_download.collect() 

    # Write results to files
    timing_file = open("timing.txt", "w")
    timing_file.write("Mapping " + str(time_mapping) + "\n")
    timing_file.write("Mapping+Downloading " + str(time_download) + "\n")
    timing_file.write("Filtering " + str(time_filtering) + "\n")
    timing_file.close()

    kmer_file = open("kmers.txt", "w")
    for item in kmers.collect():
        print>>kmer_file, item
    kmer_file.close()

    pos_file = open("positions.txt", "w")
    for item in positions.collect():
        print>>pos_file, item
    pos_file.close()

    for obj in [open(f, "r") for f in ["kmers.txt", "positions.txt"]]:
        swiftclient.client.put_object(url=storage_url, token=auth_token, container="result", name="group_11/" + str(start_time_filtering) + "_" + obj.name, contents=obj)
        obj.close()

    return

if __name__ == '__main__':
    main()
