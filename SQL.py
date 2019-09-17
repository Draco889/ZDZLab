import mysql.connector
import queue
import threading

class DATABUILD:

    def __init__(self):
        self.pipeline = queue.Queue(maxsize=10000)
        self.PATH = "/Users/Mark/PycharmProjects/ZDZLab1/Reference/"

        self.record_format = ("INSERT INTO position "
                      "(rsid, chromosome, bp) "
                      "VALUES (%s, %s, %s)")

        self.mydb = None
        self.mycursor = None

        # self.file_list = ["homo_sapiens-chr1.vcf",
        #                   "homo_sapiens-chr2.vcf",
        #                   "homo_sapiens-chr3.vcf",
        #                   "homo_sa-piens-chr4.vcf",
        #                   "homo_sapiens-chr5.vcf",
        #                   "homo_sapiens-chr6.vcf",
        #                   "homo_sapiens-chr7.vcf",
        #                   "homo_sapiens-chr8.vcf",
        #                   "homo_sapiens-chr9.vcf",
        #                   "homo_sapiens-chr10.vcf",
        #                   "homo_sapiens-chr11.vcf",
        #                   "homo_sapiens-chr12.vcf",
        #                   "homo_sapiens-chr13.vcf",
        #                   "homo_sapiens-chr14.vcf",
        #                   "homo_sapiens-chr15.vcf",
        #                   "homo_sapiens-chr16.vcf",
        #                   "homo_sapiens-chr17.vcf",
        #                   "homo_sapiens-chr18.vcf",
        #                   "homo_sapiens-chr19.vcf",
        #                   "homo_sapiens-chr20.vcf",
        #                   "homo_sapiens-chr21.vcf",
        #                   "homo_sapiens-chr22.vcf",
        #                   "homo_sapiens-chrX.vcf"]
        self.file_list  = ["homo_sapiens-chrX.vcf"]
        self.o = None
        self.log = None


    def connect(self):
        self.mydb = mysql.connector.connect(
            host="127.0.0.1",
            user="Mark",
            passwd="Draco889Tiger131",
            database="snp3"
        )

        self.mycursor = self.mydb.cursor()

    def drop_table(self):

        self.mycursor.execute("DROP TABLE position")

    def add_table(self):

        TABLE = (
            "CREATE TABLE `position` ("
            "  `rsid` bigint NOT NULL,"
            "  `chromosome` tinyint NOT NULL,"
            "  `bp` int NOT NULL"
            ") ENGINE=InnoDB")

        self.mycursor.execute(TABLE)

    def build_queue(self):
        for line in self.o:
            if "#" in line:
                continue
            else:
                split_line = line.split()
                if "rs" in split_line[2]:
                    temp = split_line[2].split("s")
                    try:
                        rsid = int(temp[1])
                        q_record = (rsid, 23, int(split_line[1]))
                        self.pipeline.put(q_record)
                    except ValueError:
                        self.log.write(line + "\n")
                else:
                    self.log.write(line + "\n")
        self.pipeline.put("stop")


    def add_record(self):
        flag = False
        count = 0
        while not flag:
            a_record = self.pipeline.get()
            if a_record != "stop":
                try:
                    self.mycursor.execute(self.record_format, a_record)
                    count += 1
                    if (count % 10000) == 0:
                        self.mydb.commit()
                        self.log.write("commit no. " + str(count % 10000) + " successful\n")
                        print("inserted " + str(count) + " records")
                except mysql.connector.errors.DataError:
                    self.log.write(str(a_record))
            else:
                self.mydb.commit()
                flag = True

    def run(self, new_table, new_log):
        self.connect()
        if new_log is True:
            self.log = open("sqllog.txt", "w+")
        else:
            self.log = open("sqllog.txt", "a")

        if new_table is True:
            self.drop_table()
            self.add_table()
        for item in self.file_list:
            self.o = open(self.PATH + item, "r")
            p = threading.Thread(target=self.build_queue)
            c = threading.Thread(target=self.add_record)
            p.start()
            c.start()

            p.join()
            c.join()
            self.o.close()
            self.log.write(item + " completed\n")
            self.log.close()
            self.log = open("sqllog.txt", "a")

        self.mycursor.close()
        self.mydb.close()

if __name__ == "__main__":

    d = DATABUILD()
    d.run(False, False)


