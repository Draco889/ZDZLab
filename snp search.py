from biomart import BiomartServer
import requests
import linecache


class snp_Search():

    def __init__(self, file, build_file, log, split_point, header_offset, save=None):
        self.search_done = False

        if save is not None:
            self.count_store = (int(linecache.getline(save, 1)) * 500 + 1)
            self.block_count = (int(linecache.getline(save, 1)) +1)
            self.log = open(log, "a")
        else:
            o = open(build_file, "w+")
            o.close()
            self.log = open(log, "w+")
            self.log.close()
            self.log = open(log, "a")
            self.count_store = header_offset
            self.block_count = 1


        while not self.search_done:
            self.server_connect()
            self.query_build(file, split_point)
            self.search(log)
            self.file_build(build_file)
            s = open("save.txt", "w+")
            s.write(str(self.block_count))
            s.close()
            self.block_count += 1

        self.results_check(file, build_file, header_offset)
        self.log.close()


    def server_connect(self):
        connected = False

        while not connected:
            try:
                server = BiomartServer("http://grch37.ensembl.org/biomart")
                connected = True
            except requests.exceptions.ConnectionError:
                connected = False
                print('Connection error')

        server.verbose = True
        self.hsapiens_snp = server.datasets['hsapiens_snp']



    def line_count(self, file, header_offset):
        self.line_count_int = 0
        flag = False
        count = 1

        while not flag:
            line = linecache.getline(file, count)
            if line == "":
                self.line_count_int = (count - header_offset)
                break
            count += 1

        return self.line_count_int


    def query_build(self, file, split_point):
        self.query_list = []
        build_done = False
        query_count = self.count_store  #first line in file is number "1"
        while not build_done:
            if query_count != (self.count_store + 500):
                line = linecache.getline(file, query_count)
                temp_list = []
                temp_string = ""
                if line != '':
                    for c in line:
                        if c is not split_point:
                            temp_list.append(c)
                        else:
                            break
                else:
                    self.search_done = True
                    break
                return_val = temp_string.join(temp_list)
                if return_val in self.query_list:
                    self.log.write("Duplicate found in block " + str(self.block_count + 1) + ": " + return_val + "\n")
                self.query_list.append(return_val)
                query_count +=1
            else:
                self.count_store = query_count
                build_done = True





    def search(self, log):
        self.data_list = []
        response = self.hsapiens_snp.search({
        'filters': {
            "snp_filter" : self.query_list
        }
        })

        for line in response.iter_lines():
            try:
                line = line.decode('utf-8')
                a = (line.split("\t"))
                self.data_list.append(str(a))

            except requests.exceptions.ConnectionError:
                continue

            except TimeoutError:
                self.log.write("Unexpected error after line " + str(self.count_store) + ", block " + str(self.block_count))
                self.log.close()

        if len(self.data_list) != len(self.query_list):
            self.log.write("Results discrepancy of " + str(len(self.query_list) - len(self.data_list)) + " in block " +
                  str(self.block_count) + "\n")

            self.log.close()
            self.log = open(log, "a")



    def file_build(self, build_file):
        self.data_count = 0
        o = open(build_file, "a")

        for item in self.data_list:
            o.write(item + "\n")
            self.data_count += 1

        o.close()


    def results_check(self, file, build_file, header_offset):
        a = self.line_count(file, header_offset)
        b = self.line_count(build_file, 0)

        if a == b:
            self.log.write("Checks Out")
        else:
            self.log.write("Count Discrepancy: " + str(a) + " (original) vs. " + str(b) + " (new)" )



s = snp_Search("GUGC_MetaAnalysis_Results_UA.csv", "test.txt", "testlog.txt", ",", 1)