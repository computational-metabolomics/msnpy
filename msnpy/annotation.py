
import os
import collections
import sqlite3
import ConfigParser
import mysql.connector
from urlparse import urlparse
import requests
import pandas as pd
from processing import mz_tolerance


class DbMolecularFormulaeLocal:

    def __init__(self, config_db):

        # Read config file
        config = ConfigParser.ConfigParser()
        config.read(config_db)

        # MySQL configurations
        # conda install -c anaconda mysql-connector-python
        self.connection = mysql.connector.connect(host=config.get('database', 'host'),
                                             passwd=config.get('database', 'password'),
                                             db=config.get('database', 'db'),
                                             user=config.get('database', 'user'),
                                             port=int(config.get('database', 'port')))
        self.cursor = self.connection.cursor()

    def table_names(self):
        ranges = []

        #try:
        #    self.cursor.execute('SELECT name FROM sqlite_master WHERE type = "table"')
        #except:
        self.cursor.execute("show tables;")

        tables = self.cursor.fetchall()
        for table in tables:
            prefixes = table[0].split("__")
            ranges.append([float(prefixes[1].replace("_", ".")), float(prefixes[2].replace("_", "."))])
        return ranges


    def select_mf(self, min_tol, max_tol, adducts=None, rules=True, sql_filter=""):
        # print min_tol, max_tol, adducts, rules
        # print (137.0588-min_tol)/(137.0588*0.000001)
        # print (137.0588-max_tol)/(137.0588*0.000001)
        # raw_input()

        if adducts is None:  # NEUTRAL LOSSES
            adducts = [None]
            adducts_lib = {None: 0.0}
        else:
            e = 0.0005486
            adducts_lib = {"[M+H]+": 1.007825 - e,
                "[M+NH4]+": 18.034374 - e,
                "[M+Na]+": 22.989770 - e,
                "[M+(39K)]+": 38.963708 - e,
                "[M+(41K)]+": 40.961825 - e,
                "[M+(6Li)]+": 6.015123 - e,
                "[M+(7Li)]+": 7.016005 - e,
                "[NaCl][M+H]+": (1.007825 + 22.989770 + 34.968853) - e,
                "[M+H+HCOOH]+": 47.01329458 - e,
                "[M+H+NaCl]+": 58.96644546 - e,
                "[M+NH4+HCOOH]+": 64.03984158 - e,
                "[M+Na+HCOOH]+": 68.99523658 - e,
                "[M+H+CHOONa]+": 68.99526858 - e,
                "[M+H+KCl]+": 74.94038516 - e,
                "[M+K+HCOOH]+": 84.96917658 - e,
                "[M-H]-": -(1.007825 - e),
                "[M+(35Cl)]-": 34.968853 + e,
                "[M+(37Cl)]-": 36.965903 + e,
                "[M+Na-2H]-": (22.989770 - (2 * 1.007825)) + e,
                "[M+K-2H]-": (38.963708 - (2 * 1.007825)) + e,
                "[M+Hac-H]-": 59.0138536}

        mf_out = []
        mf_id = 1
        table_names = self.table_names()
        for adduct in adducts:

            min_tol_temp, max_tol_temp = min_tol - adducts_lib[adduct], max_tol - adducts_lib[adduct]
            if min_tol_temp >= 0 and max_tol_temp >= 0 and min_tol_temp <= max_tol_temp:
                tables = []

                for start_end in table_names:
                    table = "mf__{}__{}".format(str(start_end[0]).replace(".", "_"), str(start_end[1]).replace(".", "_"))
                    if min_tol_temp >= start_end[0] and min_tol_temp <= start_end[1] and table not in tables:
                        tables.append(table)
                    if max_tol_temp >= start_end[0] and max_tol_temp <= start_end[1] and table not in tables:
                        tables.append(table)

                if len(tables) > 0:
                    names = ["C", "H", "N", "O", "P", "S"]
                    temp = []
                    for table in tables:
                        if rules:
                            temp.append("""select * from {} where lewis = 1 and senior = 1 and HC = 1 and NOPSC = 1 and ExactMass >= {}
                                        and ExactMass <= {} {}""".format(table, min_tol_temp, max_tol_temp, sql_filter))
                        else:
                            temp.append("select * from {} where ExactMass >= {} and ExactMass <= {}".format(table, min_tol_temp, max_tol_temp))

                    # print " UNION ".join(map(str,temp))
                    self.cursor.execute(" UNION ".join(map(str, temp)))  # + " order by ExactMass")
                    mf = self.cursor.fetchall()
                    for record in mf:
                        # print {"DBE":record[6], "LEWIS":record[7], "SENIOR":record[8], "HC":record[9], "NOPSC":record[10]}
                        mf_out.append(
                            {"mass": float(record[-1]) + adducts_lib[adduct], "atoms": dict(zip(names, record)), "adduct": adduct, "DBE": record[6],
                             "LEWIS": record[7], "SENIOR": record[8], "HC": record[9], "NOPSC": record[10]})
                        mf_id += 1
            else:
                pass
                # print "Incorrect boundaries - min: {}  max: {}".format(min_tol_temp, max_tol_temp)
                # sys.exit()

        # print min_tol, max_tol, len(mf_out)
        # print "----------------------------"

        # print([len(mf_out), min_tol, max_tol, rules])
        return mf_out


class DbMolecularFormulaeApi:

    def __init__(self, url):
        self.url = '{}/mass_range'.format(url)
        url_test = '{}//mass?mass=180.06339&tol=0.0&unit=ppm&rules=1'.format(url)
        o = urlparse(url)
        if o.scheme != "http" and o.netloc != "multiomics-int.cs.bham.ac.uk":
            raise ValueError("No database or local db available")
        else:
            r = requests.get(url_test)
            r.raise_for_status()

    def select_mf(self, min_tol, max_tol, adducts=None, rules=True, sql_filter=""):
        # print min_tol, max_tol, adducts, rules
        # print (137.0588-min_tol)/(137.0588*0.000001)
        # print (137.0588-max_tol)/(137.0588*0.000001)
        # raw_input()

        if adducts is None:  # NEUTRAL LOSSES
            adducts = [None]
            adducts_lib = {None: 0.0}
        else:
            e = 0.0005486
            adducts_lib = {"[M+H]+": 1.007825 - e,
                "[M+NH4]+": 18.034374 - e,
                "[M+Na]+": 22.989770 - e,
                "[M+(39K)]+": 38.963708 - e,
                "[M+(41K)]+": 40.961825 - e,
                "[M+(6Li)]+": 6.015123 - e,
                "[M+(7Li)]+": 7.016005 - e,
                "[NaCl][M+H]+": (1.007825 + 22.989770 + 34.968853) - e,
                "[M+H+HCOOH]+": 47.01329458 - e,
                "[M+H+NaCl]+": 58.96644546 - e,
                "[M+NH4+HCOOH]+": 64.03984158 - e,
                "[M+Na+HCOOH]+": 68.99523658 - e,
                "[M+H+CHOONa]+": 68.99526858 - e,
                "[M+H+KCl]+": 74.94038516 - e,
                "[M+K+HCOOH]+": 84.96917658 - e,
                "[M-H]-": -(1.007825 - e),
                "[M+(35Cl)]-": 34.968853 + e,
                "[M+(37Cl)]-": 36.965903 + e,
                "[M+Na-2H]-": (22.989770 - (2 * 1.007825)) + e,
                "[M+K-2H]-": (38.963708 - (2 * 1.007825)) + e,
                "[M+Hac-H]-": 59.0138536}

        mf_out = []
        mf_id = 1
        for adduct in adducts:

            params = {"lower": min_tol - adducts_lib[adduct],
                      "upper": max_tol - adducts_lib[adduct], "rules": int(rules)}
            response = requests.get(self.url, params=params)

            # Check response is ok
            if not response:
                continue
            resp_d = response.json()

            # check records key in response dict/json
            if not 'records' in resp_d:
                continue
            records = resp_d["records"]

            names = ["C", "H", "N", "O", "P", "S"]
            for record in records:
                
                # print {"DBE":record[6], "LEWIS":record[7], "SENIOR":record[8], "HC":record[9], "NOPSC":record[10]}
                mf_out.append(
                    {"mass": float(record['ExactMass']) + adducts_lib[adduct],
                     "atoms": {k: record[k] for k in names},
                     "adduct": adduct,
                     "DBE": record['DoubleBondEquivalents'],
                     "LEWIS": record['LEWIS'],
                     "SENIOR": record['SENIOR'],
                     "HC": record['HC'],
                     "NOPSC": record['NOPSC']})
                mf_id += 1
            else:
                pass
                # print "Incorrect boundaries - min: {}  max: {}".format(min_tol_temp, max_tol_temp)
                # sys.exit()

        # print min_tol, max_tol, len(mf_out)
        # print "----------------------------"

        # print([len(mf_out), min_tol, max_tol, rules])
        return mf_out


def mz_pair_diff_tol(lower_mz, upper_mz, lower_mz_tol, upper_mz_tol, lower_mz_unit="ppm", upper_mz_unit="ppm"):
    mz_diff = upper_mz - lower_mz
    lmt = mz_tolerance(lower_mz, lower_mz_tol, lower_mz_unit)[1] - lower_mz
    hmt = mz_tolerance(upper_mz, upper_mz_tol, upper_mz_unit)[1] - upper_mz
    total_tol = lmt + hmt

    if mz_diff - total_tol < 0.0:
        return 0.0, mz_diff + total_tol
    else:
        return mz_diff - total_tol, mz_diff + total_tol


def annotate_mf(spectral_trees, db_out, ppm, adducts=["[M+H]+"], rules=True, mf_db="http://multiomics-int.cs.bham.ac.uk", prefix_inp=""):

    if os.path.isfile(mf_db):
        db = DbMolecularFormulaeLocal(config_db=mf_db)
    else:
        db = DbMolecularFormulaeApi(url=mf_db)

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    for G in spectral_trees:
        if prefix_inp == "":
            prefix = "_{}".format(G.graph["id"])
        else:
            prefix = "_{}_{}".format(G.graph["id"], prefix_inp)

        cursor.execute("""DROP TABLE IF EXISTS MF{}""".format(prefix))
        cursor.execute("""
        CREATE TABLE MF{} (
            MZ_ID       TEXT,
            MF_ID       INTEGER,
            C           INTEGER,
            H           INTEGER,
            N           INTEGER,
            O           INTEGER,
            P           INTEGER,
            S           INTEGER,
            DBE         INTEGER,
            LEWIS       INTEGER,
            SENIOR      INTEGER,
            HC          INTEGER,
            NOPSC       INTEGER,
            ADDUCT      TEXT,
            MASS        FLOAT,
            MZ          FLOAT,
            PPM_ERROR   FLOAT,
            PRECURSOR   INTEGER,
            MSLEVEL     TEXT,
            FLAG        INTEGER,
            PRIMARY KEY (MZ_ID, MF_ID)
        )""".format(prefix))

        mf_id = 1
        nodes_done = []
        rows = []
        for edge in G.edges(data=True):

            if edge[2]['mzdiff'] > 0.5 and edge[2]['type'] == "e":

                for node_id in [edge[0], edge[1]]:

                    node = G.node[node_id]

                    if node_id not in nodes_done:

                        nodes_done.append(node_id)
                        min_tol, max_tol = mz_tolerance(node["mz"], ppm)  # ppm tol

                        records_mf = db.select_mf(min_tol, max_tol, adducts, rules)

                        if len(records_mf) > 0:
                            for mf in records_mf:
                                ppm_error = round((node["mz"] - mf["mass"]) / (mf["mass"] * 0.000001), 2)
                                values = (
                                node_id, mf_id, mf["atoms"]["C"], mf["atoms"]["H"], mf["atoms"]["N"], mf["atoms"]["O"],
                                mf["atoms"]["P"], mf["atoms"]["S"], mf["adduct"], mf["mass"], node["mz"], ppm_error, int(node["precursor"]),
                                node["mslevel"], mf["DBE"], mf["LEWIS"], mf["SENIOR"], mf["HC"], mf["NOPSC"], 0,)
                                rows.append(values)
                                mf_id += 1
                        else:
                            values = (node_id, None, None, None, None, None, None, None, None, None, node["mz"], None,
                                      int(node["precursor"]), node["mslevel"], None, None, None, None, None, None,)
                            rows.append(values)

                node_i = G.node[edge[0]]
                node_j = G.node[edge[1]]

                min_tol, max_tol = mz_pair_diff_tol(node_j["mz"], node_i["mz"], ppm, ppm)

                records_mf = db.select_mf(min_tol, max_tol, None, False)

                if len(records_mf) > 0:
                    for mf in records_mf:
                        ppm_error = round((mf["mass"] - edge[2]['mzdiff']) / (mf["mass"] * 0.000001), 2)
                        values = ("{}__{}".format(edge[0], edge[1]), mf_id,
                                  mf["atoms"]["C"], mf["atoms"]["H"],
                                  mf["atoms"]["N"], mf["atoms"]["O"],
                                  mf["atoms"]["P"], mf["atoms"]["S"],
                                  None, mf["mass"], edge[2]['mzdiff'], ppm_error,
                                  None,
                                  "{}_{}".format(node_i["mslevel"], node_j["mslevel"]),
                                  mf["DBE"], mf["LEWIS"], mf["SENIOR"], mf["HC"], mf["NOPSC"], 0)
                        rows.append(values)
                        mf_id += 1
                else:
                    values = (
                    "{}__{}".format(edge[0], edge[1]), None, None, None, None, None, None, None, None, None, None,
                    edge[2]['mzdiff'], None, "{}_{}".format(node_i["mslevel"], node_j["mslevel"]),
                    None, None, None, None, None, 0,)
                    rows.append(values)
                    ####################################
                    # else:
                    #    values = ("{}__{}".format(edge[0], edge[1]), None, None, None, None, None, None, None, None, None, None,
                    #              edge[2]['mzdiff'], None, "{}_{}".format(node_i["mslevel"], node_j["mslevel"]),
                    #              None, None, None, None, None, 0,)
                    #    rows.append(values)

        cursor.executemany("""
        INSERT INTO MF{}
               (MZ_ID, MF_ID, C, H, N, O, P, S, ADDUCT, MASS, MZ, PPM_ERROR, PRECURSOR, MSLEVEL, DBE, LEWIS, SENIOR, HC, NOPSC, flag)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """.format(prefix), rows)

        cursor.execute("""
        CREATE INDEX CHNOPS_ADDUCT_{}_IDX
            ON MF{} (MZ_ID, PRECURSOR, ADDUCT, C, H, N, O, P, S, FLAG);""".format(prefix, prefix))
        cursor.execute("""
        CREATE INDEX MZ_ID_{}_IDX
            ON MF{} (MZ_ID);""".format(prefix, prefix))
        cursor.execute("""
        CREATE INDEX MF_MSL_P_F_{}_IDX
            ON MF{} (MSLEVEL, PRECURSOR, FLAG)""".format(prefix, prefix))
        cursor.execute("""
        CREATE INDEX MF_FLAG_{}_IDX
            ON MF{} (FLAG)""".format(prefix, prefix))
        conn.commit()

    return spectral_trees


def mf_tree(G, path_db, max_mslevel, prefix):

    conn = sqlite3.connect(path_db)
    cursor = conn.cursor()

    mft = []

    cursor.execute("""
    SELECT MF_ID
      FROM MF{}
     WHERE FLAG > 0 AND PRECURSOR = 1 AND MSLEVEL = 1
    """.format(prefix))
    mfs_precs = cursor.fetchall()

    for mf_prec in mfs_precs:

        sql_str = "SELECT DISTINCT \n".format(prefix)
        sql_str += "E1.MZ_ID_PREC, E1.MZ_ID_NL, E1.MZ_ID_FRAG, E1.MF_ID_PREC, E1.MF_ID_NL, E1.MF_ID_FRAG,\n"
        sql_col = []
        sql_join = []

        for i in range(2, max_mslevel + 1):
            sql_col.append("E{}.MZ_ID_PREC, E{}.MZ_ID_NL, E{}.MZ_ID_FRAG, E{}.MF_ID_PREC, E{}.MF_ID_NL, E{}.MF_ID_FRAG".format(i, i, i, i, i, i))
            sql_join.append("LEFT OUTER JOIN\nedges{} AS E{}\nON E{}.MF_ID_FRAG = E{}.MF_ID_PREC".format(prefix, i, i - 1, i))

        sql_str += ",\n".join(map(str, sql_col))
        sql_str += "\nFROM edges{} AS E1\n".format(prefix)
        sql_str += "\n".join(map(str, sql_join))
        sql_str += "\nWHERE E1.MF_ID_PREC = {};".format(mf_prec[0]) # ID

        cursor.execute(sql_str)
        records = cursor.fetchall()

        GG = G.copy()
        atoms = ["C", "H", "N", "O", "P", "S"]
        for record in records:
            record = filter(lambda x: x != None, record)
            for i in range(0, len(record), 6):
                ids = record[i:i + 6]
                cursor.execute("""
                SELECT MZ_ID, MF_ID, C, H, N, O, P, S, ADDUCT, MASS
                  FROM MF{}
                 WHERE MF_ID = {} OR MF_ID = {} and FLAG > 0
                """.format(prefix, ids[3], ids[5]))

                for mf in cursor.fetchall():
                    if "mf" not in GG.node[str(mf[0])]:
                        GG.node[str(mf[0])]["mf"] = {}
                    if mf[1] not in GG.node[str(mf[0])]["mf"]:
                        #GG.node[str(mf[0])]["mf"][str(mf[1])] = {"mass": float(mf[14]), "atoms": collections.OrderedDict(zip(atoms, mf[2:8])), "ion": mf[13]}
                        GG.node[str(mf[0])]["mf"][str(mf[1])] = {"mass": float(mf[9]), "mf": print_formula(collections.OrderedDict(zip(atoms, mf[2:8]))), "adduct": mf[8]}

                cursor.execute("""
                SELECT MZ_ID, MF_ID, C, H, N, O, P, S, ADDUCT, MASS
                  FROM MF{}
                 WHERE MF_ID = '{}' and FLAG > 0
                """.format(prefix, ids[4]))

                for mf in cursor.fetchall():
                    n = mf[0].split("__")
                    if "mf" not in GG[n[0]][n[1]]:
                        GG[n[0]][n[1]]["mf"] = {}
                    if mf[1] not in GG[n[0]][n[1]]["mf"]:
                        GG[n[0]][n[1]]["mf"][str(mf[1])] = {"mass": float(mf[9]), "mf": print_formula(collections.OrderedDict(zip(atoms, mf[2:8])))}

        #print "N:", GG.number_of_nodes(), "E:", GG.number_of_edges()
        nodes = GG.nodes(data=True)
        for node in nodes:
            if "mf" not in node[1]:
                GG.remove_node(node[0])
        GG.graph["id"] = "{}_{}".format(GG.graph["id"], mf_prec[0])
        mft.append(GG)

    conn.close()
    return mft


def filter_mf(trees, path_db):

    #http://www.sqlstyle.guide/

    conn = sqlite3.connect(path_db)
    cursor = conn.cursor()

    atoms = ["C", "H", "N", "O", "P", "S"]

    annotated_trees = []
    for G in trees:

        prefix = "_{}".format(G.graph["id"])

        ################################################################
        # TABLES: EDGES AND MZ_PREC_FRAG
        ################################################################
        cursor.execute("""DROP TABLE IF EXISTS EDGES{}""".format(prefix))
        cursor.execute("""
        CREATE TABLE EDGES{} (
            MF_ID_PREC      INTERGER,
            MF_ID_NL        INTERGER,
            MF_ID_FRAG      INTERGER,
            MF_ID_PREC_MASS FLOAT,
            MF_ID_NL_MASS   FLOAT,
            MF_ID_FRAG_MASS FLOAT,
            MZ_ID_PREC      TEXT,
            MZ_ID_NL        TEXT,
            MZ_ID_FRAG      TEXT,
            MF_ID_PREC_FLAG INTEGER,
            MF_ID_NL_FLAG   INTEGER,
            MF_ID_FRAG_FLAG INTERGER,
            PRIMARY KEY(MF_ID_PREC, MF_ID_NL, MF_ID_FRAG)
        )""".format(prefix))

        cursor.execute("""
        CREATE INDEX EDGES_IDS{}_IDX
            ON EDGES{} (MF_ID_PREC, MF_ID_NL, MF_ID_FRAG);
        """.format(prefix, prefix))

        cursor.execute("""
        CREATE INDEX EDGES_FLAGS{}_IDX
            ON EDGES{} (MF_ID_PREC_FLAG, MF_ID_NL_FLAG, MF_ID_FRAG_FLAG)
        """.format(prefix, prefix))

        cursor.execute("""DROP TABLE IF EXISTS MZ_PREC_FRAG{}""".format(prefix))

        cursor.execute("""
        CREATE TABLE MZ_PREC_FRAG{} (
            MZ_ID_PREC   INTERGER,
            MZ_ID_FRAG   INTERGER,
            PRIMARY KEY (MZ_ID_PREC, MZ_ID_FRAG)
        )""".format(prefix))
        conn.commit()

        cursor.execute("""
        SELECT COUNT(MF_ID)
          FROM MF{}
         WHERE FLAG >= 0
        """.format(prefix))

        improvement = [cursor.fetchone()[0]]

        if G.number_of_nodes() == 0:
            continue

        # Insert all inital edges
        for edge in [edge for edge in G.edges(data=True) if edge[2]["type"] == "e" and edge[2]["mzdiff"] > 0.5]:

            # Set IDs to search sqlite database
            mz_id_prec, mz_id_frag, mz_id_nl = edge[0], edge[1], "{}__{}".format(edge[0], edge[1])

            sub_queries = ['MF1.MZ_ID = "{}" AND MF2.MZ_ID = "{}" AND MF3.MZ_ID = "{}" '
                           'AND MF1.PRECURSOR = 1 AND MF1.ADDUCT = MF3.ADDUCT'.format(mz_id_prec, mz_id_nl, mz_id_frag),
                           'MF1.FLAG >= 0 AND MF2.FLAG >= 0 AND MF3.FLAG >= 0']

            # Number of C in the fragment MF should be smaller or equal to the number of C in the precursor MF
            # Number of C in the fragment MF plus the number of C in the neutral mass should be equal to the number of C in the precursor MF
            # Apply for N, O, P and S.
            for atom in atoms:
                sub_queries.append("MF3.{} <= MF1.{} AND MF3.{} + MF2.{} = MF1.{}".format(atom, atom, atom, atom, atom))

            # APPLY CONSTRAINS & INSERT EDGES
            cursor.execute("""
            INSERT INTO MZ_PREC_FRAG{} (MZ_ID_PREC, MZ_ID_FRAG)
            VALUES ('{}','{}')""".format(prefix, mz_id_prec, mz_id_frag))

            cursor.execute("""
            INSERT INTO EDGES{}
                   SELECT MF1.MF_ID, MF2.MF_ID, MF3.MF_ID, MF1.MASS, MF2.MASS, MF3.MASS, MF1.MZ_ID, MF2.MZ_ID, MF3.MZ_ID, 0, 0, 0
                     FROM MF{} as MF1, MF{} as MF2, MF{} as MF3
                    WHERE {}
            """.format(prefix, prefix, prefix, prefix," AND ".join(map(str, sub_queries))))
        conn.commit()

        cursor.execute("""
        SELECT MAX(CAST(mslevel as int))
          FROM MF{}
         WHERE mslevel
           NOT LIKE '%\_%'
        """.format(prefix))
        max_mslevel = cursor.fetchone()[0]

        if max_mslevel == None: max_mslevel = 1

        for loop in range(1, max_mslevel + 1):

            # loop through neutral losses (i.e. type is e)
            for edge in [edge for edge in G.edges(data=True) if edge[2]["type"] == "e" and edge[2]["mzdiff"] > 0.5]:

                # Define constrains as above
                mz_id_prec, mz_id_frag, mz_id_nl = edge[0], edge[1], "{}__{}".format(edge[0], edge[1])

                sub_queries = ['MF1.MZ_ID = "{}" AND MF2.MZ_ID = "{}" AND MF3.MZ_ID = "{}" '
                    'AND MF1.PRECURSOR = 1 AND MF1.ADDUCT = MF3.ADDUCT'.format(mz_id_prec, mz_id_nl, mz_id_frag),
                    'MF1.FLAG >= 0 AND MF2.FLAG >= 0 AND MF3.FLAG >= 0']

                for atom in atoms:
                    sub_queries.append("MF3.{} <= MF1.{} AND MF3.{} + MF2.{} = MF1.{}".format(atom, atom, atom, atom, atom))

                # Apply constrains & Update edges
                cursor.execute("""
                UPDATE EDGES{}
                   SET MF_ID_PREC_FLAG = {}, MF_ID_NL_FLAG = {}, MF_ID_FRAG_FLAG = {}
                 WHERE
                    EXISTS(
                        SELECT MF1.MF_ID, MF2.MF_ID, MF3.MF_ID
                          FROM MF{} AS MF1, MF{} as MF2, MF{} AS MF3
                         WHERE {}
                           AND EDGES{}.MF_ID_PREC = MF1.MF_ID
                           AND EDGES{}.MF_ID_NL = MF2.MF_ID
                           AND EDGES{}.MF_ID_FRAG = MF3.MF_ID
                    )
                """.format(prefix, loop, loop, loop, prefix, prefix, prefix, " AND ".join(map(str, sub_queries)), prefix, prefix, prefix))  # AND MF1.FLAG >= 1 AND MF1.FLAG > 1  AND MF1.FLAG >= 1;
            conn.commit()

            conn.commit()

            # Update flag MF table based on the flag in table edges and the number of loops
            values = (loop, loop, loop, loop,)
            cursor.execute("""
            UPDATE MF{}
               SET FLAG = ?
             WHERE MF_ID 
                IN (
                   SELECT MF_ID_PREC
                     FROM EDGES{}
                    WHERE MF_ID_PREC_FLAG = ?
                )
                OR MF_ID 
                IN (
                   SELECT MF_ID_NL
                     FROM EDGES{}
                    WHERE MF_ID_NL_FLAG = ?
                )
                OR MF_ID 
                IN (
                    SELECT MF_ID_FRAG
                      FROM EDGES{}
                     WHERE MF_ID_PREC_FLAG = ?
                )
            """.format(prefix, prefix, prefix, prefix), values)

            ###########################################################################
            # 1a. Select all fragment MF that do not occur as a valid precursor or are part of a 'loose' edge
            # Example: C9H6NO Not valid for all MF rules, see http://link.springer.com/article/10.1007/s11419-016-0319-8
            # Example: C9H6NO Not valid for all MF rules, see http://metabolomics.jp/wiki/Ojima:KOX00444n

            # pos_5-methoxy-3-indoleacetic_acid__run_1

            # No annotation for lower MS level fragments
            # pos_DL_3_indolelactic_acid__run_1
            # e.g.m/z 170.06 > 142.07 > 115.05 > 67.70 (electrical noise peak)

            # Count number of edges related to a fragment
            # Should have a minimum of two edges if not a MS level 1 precursor

            cursor.execute("""
            SELECT MF_ID, MZ_ID, count(*) as C
              FROM (
                    SELECT DISTINCT E1.MF_ID_PREC AS MF_ID, E1.MZ_ID_PREC AS MZ_ID
                      FROM EDGES{} AS E1
                     WHERE E1.MF_ID_PREC_FLAG = ? AND E1.MF_ID_FRAG_FLAG = ?
                     UNION ALL
                    SELECT DISTINCT E2.MF_ID_FRAG AS MF_ID, E2.MZ_ID_FRAG AS MZ_ID
                      FROM EDGES{} AS E2
                     WHERE E2.MF_ID_PREC_FLAG = ? AND E2.MF_ID_FRAG_FLAG = ?
                   )
             WHERE MF_ID NOT 
                IN (
                    SELECT MF_ID
                      FROM MF{}
                     WHERE MSLEVEL = 1
                   )
               AND MZ_ID 
                IN (
                    SELECT MZ_ID_PREC
                      FROM MZ_PREC_FRAG{}
                   )
             GROUP BY MF_ID
            HAVING C < 2
            """.format(prefix, prefix, prefix, prefix), (loop, loop, loop, loop, ))
            records = cursor.fetchall()
            for record in records:

                # Set the flag for the fragment/precursor MF to a minus value (i.e. -flag)
                # If all related fragment masses (if exists) have no MF annotation

                values = (-loop, record[0], loop, )
                cursor.execute("""
                UPDATE MF{}
                   SET FLAG = ?
                 WHERE MF{}.MF_ID = ? AND MF{}.PRECURSOR = 1 AND MF{}.FLAG = ?
                """.format(prefix, prefix, prefix, prefix, prefix), values)

            conn.commit()
            # UPDATE COLUMN FLAG IN TABLE EDGES FOR PERCURSOR, NEUTRAL LOSS AND FRAGMENT BASED ON MF TABLE

            cursor.execute("""
            UPDATE EDGES{}
               SET MF_ID_FRAG_FLAG = (SELECT FLAG FROM MF{} WHERE MF_ID = MF_ID_FRAG)
             WHERE
            EXISTS (
                    SELECT FLAG 
                      FROM MF{}
                     WHERE MF_ID = MF_ID_FRAG
                   );
            """.format(prefix, prefix, prefix))

            cursor.execute("""
            UPDATE EDGES{}
               SET MF_ID_NL_FLAG = (
                                    SELECT FLAG 
                                    FROM MF{}
                                    WHERE MF_ID = MF_ID_NL
                                   )
             WHERE
            EXISTS (
                    SELECT FLAG 
                      FROM MF{}
                     WHERE MF_ID = MF_ID_NL
                   );
            """.format(prefix, prefix, prefix))

            cursor.execute("""
            UPDATE EDGES{}
               SET MF_ID_PREC_FLAG = (
                                      SELECT FLAG 
                                      FROM MF{}
                                      WHERE MF_ID = MF_ID_PREC
                                     )
             WHERE
            EXISTS(
                   SELECT FLAG
                     FROM MF{}
                    WHERE MF_ID = MF_ID_PREC
                  );
            """.format(prefix, prefix, prefix))
            conn.commit()
            ###################################################################################################

            values = (-loop, loop - 1, )
            cursor.execute("""
            UPDATE MF{}
               SET FLAG = ?
             WHERE FLAG = ?""".format(prefix), values)
            conn.commit()

            cursor.execute("""
            SELECT COUNT(MF_ID)
              FROM MF{}
             WHERE FLAG >= 0
            """.format(prefix))
            improvement.append(cursor.fetchone()[0])

        print "MF AFTER CONSTRAINS: ", ", ".join(map(str, improvement))
        annotated_trees.extend(mf_tree(G, path_db, max_mslevel, prefix))

    conn.close()
    return annotated_trees


def print_formula(atom_counts):
    formula_out = ""

    order_atoms = ["C", "H", "N", "O", "P", "S"]

    for atom in order_atoms:
        if atom_counts[atom] > 1:
            formula_out += atom + str(atom_counts[atom])
        elif abs(atom_counts[atom]) == 1:
            formula_out += atom
        else:
            pass

    return formula_out


def rank_mf(trees):

    columns = ["TreeID", 'GroupID', 'MolecularFormulaID', 'MolecularFormula', 'Adduct', 'Rank', 'TotalRanks', 'RankedEqual', 'Trees', 'NeutralLossesExplained']
    df = pd.DataFrame(columns=columns)

    annotated_trees = collections.OrderedDict()
    for graph in trees:
        tree_id = graph.graph["id"].split("_")
        annotated_trees.setdefault(tree_id[0], []).append(graph)

    for i, graphs in enumerate(annotated_trees.values()):

        if len(graphs) == 0:
            continue
        df_subset = pd.DataFrame(columns=columns)
        for graph in graphs:
            mf_id = graph.graph["id"].split("_")[1]
            group_id = graph.graph["id"].split("_")[0]
            mf = str(graph.node[list(graph.nodes())[0]]["mf"][str(mf_id)]["mf"])
            adduct = str(graph.node[list(graph.nodes())[0]]["mf"][str(mf_id)]["adduct"])
            values = [graph.graph["id"], group_id, mf_id, mf, adduct, 0, 0, 0, len(graphs), graph.number_of_edges()]#, mf_str, ion_str]
            d = collections.OrderedDict(zip(columns, values))
            df_subset = df_subset.append(d, ignore_index=True)

        df_subset["Rank"] = df_subset['NeutralLossesExplained'].rank(method='dense', ascending=False).astype(int)
        df_subset['RankedEqual'] = df_subset.groupby('NeutralLossesExplained')['NeutralLossesExplained'].transform('count')
        df_subset['TotalRanks'] = df_subset['Rank'].nunique()
        df_subset = df_subset.sort_values(by=['Rank', 'MolecularFormulaID'])
        df = pd.concat([df, df_subset], ignore_index=True)
    return df
