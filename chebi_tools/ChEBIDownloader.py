import os
import gzip
import shutil
import wget
import logging

import pandas as pd

from pathlib import Path as P


class ChEBIDownloader:
    def __init__(self, download_dir=None):
        """
        Initializes the class and sets the default download directory.
        If download_dir is not provided, the LIBCHEBIPY_DOWNLOAD_DIR
        environment variable or a default directory is used.
        """
        self.download_dir = (
            download_dir
            if (download_dir is not None)
            else os.environ.get(
                "LIBCHEBIPY_DOWNLOAD_DIR",
                os.path.join(os.path.expanduser("~"), "libChEBI"),
            )
        )

        self.url = "https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/"
        self.url_obo = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/"

        self.fns = {
            "chebiId_inchi": dict(
                orig="chebiId_inchi.tsv",
                processed=P(self.download_dir) / "chebiId_inchi.parquet",
            ),
            "comments": dict(
                orig="comments.tsv", processed=P(self.download_dir) / "comments.parquet"
            ),
            "chemical_data": dict(
                orig="chemical_data.tsv",
                processed=P(self.download_dir) / "chemical_data.parquet",
            ),
            "compounds": dict(
                orig="compounds.tsv.gz",
                processed=P(self.download_dir) / "compounds.parquet",
            ),
            "database_accession": dict(
                orig="database_accession.tsv",
                processed=P(self.download_dir) / "database_accession.parquet",
            ),
            "names": dict(
                orig="names.tsv.gz", processed=P(self.download_dir) / "names.parquet"
            ),
            "reference": dict(
                orig="reference.tsv.gz",
                processed=P(self.download_dir) / "reference.parquet",
            ),
            "relation": dict(
                orig="relation.tsv", processed=P(self.download_dir) / "relation.parquet"
            ),
            "structures": dict(
                orig="structures.csv.gz",
                processed=P(self.download_dir) / "structures.parquet",
            ),
            "obo": dict(
                orig="chebi_core.obo.gz",
                processed=P(self.download_dir) / "chebi_core.obo",
            ),
        }

    def download(self, what="compounds"):
        """
        Downloads a specific file from ChEBI to the download directory.
        :param what: (str) the file to download, must be one of the keys in self.fns
        """
        P(self.download_dir).mkdir(parents=True, exist_ok=True)
        if what != "obo":
            url = self.url + self.fns[what]["orig"]
        else:
            url = self.url_obo + self.fns[what]["orig"]
        print(f"Downloading {url} to {self.download_dir}")
        fn_target = P(self.download_dir) / self.fns[what]["orig"]
        wget.download(url, str(fn_target))
        if fn_target.suffix == ".gz":
            fn_target = self.extract(fn_target, fn_target.with_suffix(""))
        self.convert_to_parquet(fn_target)

    def download_all(self):
        """
        Downloads all files from ChEBI to the download directory.
        """
        for what in self.fns.keys():
            self.download(what=what)
        print("Done")

    def extract(self, fn, fn_out):
        """
        Extracts a gzip file to a new file.
        :param fn: (str or Path) the gzip file to extract
        :param fn_out: (str or Path) the file to write the extracted data to
        :return: (str or Path) the extracted file
        """
        print(f"Extracting: {fn} to {fn_out}")
        with gzip.open(fn, "rb") as f_in:
            with open(fn_out, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        if P(fn_out).is_file():
            P(fn).unlink()
        return fn_out

    def convert_to_parquet(self, fn):
        """
        Converts a file to the parquet format.
        :param fn: (str or Path) the file to convert
        :return: (str or Path) the converted file
        """
        fn = P(fn)
        fn_out = fn.with_suffix(".parquet")
        print(f"Converting: {fn} to {fn_out}")
        if fn.name == "structures.csv":
            pd.read_csv(fn, na_filter=False, low_memory=False).to_parquet(fn_out)
        elif fn.suffix == ".tsv":
            pd.read_csv(
                fn,
                sep=None if fn.suffix == "csv" else "\t",
                na_filter=False,
                encoding="437",
                low_memory=False,
            ).to_parquet(fn_out)
        else:
            return fn
        assert fn_out.is_file(), fn_out
        fn.unlink()
        return fn_out

    def check_files(self):
        """
        Check whether the file has been downloaded
        :param fn: (str or Path) the file to check
        :return: (bool) whether the file is complete or not
        """
        missing = [fn for fn, path in self.files.items() if not path.is_file()]
        return missing

    def download_missing(self):
        """
        Downloads missing files from ChEBI to the download directory.
        """
        missing = self.check_files()
        if missing:
            for what in missing:
                print(what)
                self.download(what=what)

    @property
    def files(self):
        return {fn: self.fns[fn]["processed"] for fn in self.fns}

    def get_path(self, name):
        return self.files[name]
