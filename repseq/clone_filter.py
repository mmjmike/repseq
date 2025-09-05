from .clonosets import (decide_count_and_frac_columns,
                        get_column_names_from_clonoset,
                        filter_by_functionality)

from .common_functions import (extract_segment,
                               extract_refpoint_position, overlap_type_to_flags)

from collections.abc import Hashable

import random

class Filter:
    
    """
    Clonoset filter.
    May be used to filter clonosets:
        - by clone functionality
        - randomly downsample them to particular number of reads of UMIs
        - take top clonotypes by size (number of reads of UMIs) with or without random mixing

    Args:
        name (str): the name of the filter. Will be displayed in print
        functionality (str): Possible values:
            - "a" - any (default). No clones are filtered out
            - "f" - only functional. Those, not having stop codons and 
                frame shifts in CDR3 regions, or having non-empty values in CDR3 amino-acid
                sequence
            - "n" - only-nonfunctional - opposite to "f" - functional
        downsample (int): the number of reads/UMIs to randomly downsample the clonoset to.
            default value 'None' - means not to apply downsampling
        top (int): the number of top biggest by reads/UMIs clonotypes to take from the clonoset.
            default value 'None' - means not to apply top
        by_umi (bool): default=False. Which column to take for clonotype count - reads or UMIs 
            (if UMI count column exists).
        mix_tails (bool): default=False. Defines whether to randomly mix-up the order of clonotypes
            before sorting by size and taking the top clonotypes. Basically mix_tails=True mixes up 
            clonotypes with the same size in read or UMIs.
        count_threshold (int): limits [0:100000], all clonotypes with count less than this value will
            be filtered out
        seed (any hashable type): better to use int - seed for reproducibility of random events 
            (downsampling or top with mix-tails). Default=None.
        unweight (bool): each clonotype counts (either reads or UMIs) are set to 1
        recount_fractions (bool): if `True`, clonotype fractions are recalculated after filtration
        white_list (list of tuples): If specified, only clonotypes matching those listed will be retained. Either `aa`, `aaV` or `aaVJ` 
            formats can be used to list clonotypes, e.g. [(“CASSS..”)], [(“CASSS..”, “TRBV2”)] or [(“CASSS..”, “TRBV2”, “TRBJ1”)]. 
            It is applied before `black_list`
        black_list (list of tuples): If specified, only clonotypes not listed will be retained. Either `aa`, `aaV` or `aaVJ` 
            formats can be used to list clonotypes, e.g. [(“CASSS..”)], [(“CASSS..”, “TRBV2”)] or [(“CASSS..”, “TRBV2”, “TRBJ1”)]
        pool_clonoset_by (str): possible values are ["", "aa", "aaV", "aaVj", "nt", "ntV", "ntVJ"]. Clones with identical parameters are merged, 
            keeping the largest one, while their counts are summed.
        convert (bool): By default, columns are added to the clonotype set to convert 
            it to the VDJtools format, and the original columns are removed. If set to  `False`, all original columns are preserved.
        ignore_small_clonosets (bool): If the top or downsample threshold exceeds the counts for a clonoset, the set is kept intact
    """

    def __init__(self, name="default_filter", functionality="a", downsample=None,
                 top=None, by_umi=False, mix_tails=False, count_threshold=None, 
                 unweight=False, seed=None, recount_fractions=True,
                 white_list=[], black_list=[], pool_clonoset_by="", convert=True, 
                 ignore_small_clonosets=False):
        self.name = name
        self.functionality = functionality
        self.downsample_size = downsample
        self.top = top
        self.by_umi = by_umi
        self.mix_tails = mix_tails
        self.seed = seed
        self.count_threshold = count_threshold
        self.unweight = unweight
        self.recount_fractions = recount_fractions
        self.white_list = white_list
        self.black_list = black_list
        self.pool_by = pool_clonoset_by
        self.convert = convert
        self.ignore_small_clonosets = ignore_small_clonosets
        self._check_input()
        
    def spawn(self):
        """
        
        Returns:
            the copy of the filter. Necessary for parallel computing

        """
        return Filter(name=self.name, functionality=self.functionality,
                      downsample=self.downsample_size, top=self.top,
                      by_umi=self.by_umi, mix_tails=self.mix_tails,
                      count_threshold=self.count_threshold, seed=self.seed,
                      unweight=self.unweight,
                      recount_fractions=self.recount_fractions,
                      white_list = self.white_list,
                      black_list = self.black_list,
                      ignore_small_clonosets=self.ignore_small_clonosets
                      )
        
    def apply(self, input_clonoset, colnames=None):
        """
        Main method of the Filter object - application of it to a clonoset

        Args:
            input_clonoset (pd.DataFrame): clonoset in the form of Pandas DataFrame in
                MiXCR(3 or 4+ version), VDJtools or Bioadaptive formats.
            colnames (dict, optional): Dictionary of available specific column names.
                Defaults to None - colnames imputed automatically.

        Returns:
            clonoset (pd.DataFrame): clonoset after converting to common (VDJtools-like)
                format and applying functionality filtration and downsampling or taking top
        """
        
        # copy clonoset for not changing the original one
        clonoset = input_clonoset.copy()
        if colnames is None:
            colnames = get_column_names_from_clonoset(clonoset)

        # create common columns: vdj-refPoints and VDJC-segments in common state
        clonoset = self._make_common_columns(clonoset, colnames)
        colnames = get_column_names_from_clonoset(clonoset)

        # converting to common VDJtools-like format and obtaining new colnames
        if self.convert:
            clonoset = self._convert_clonoset(clonoset, colnames)
            colnames = get_column_names_from_clonoset(clonoset)

        # application of main filters
        if self.functionality != "a":
            clonoset = self._filter_by_functionality(clonoset, colnames)
        
        clonoset = self._filter_by_count(clonoset, colnames)
        
        clonoset = self._downsample(clonoset, colnames)
        clonoset = self._get_top(clonoset, colnames)

        if self.unweight:
            clonoset = self._unweight(clonoset, colnames)
        # the fraction columns need to be recounted after filtering, as they
        # remain the same as in the original clonoset before filtration
        if self.recount_fractions:
            clonoset = self._recount_fractions_for_clonoset(clonoset, colnames)
        if self.pool_by:
            clonoset = self._pool_clonoset(clonoset, colnames)
        if len(self.white_list) > 0:
            clonoset = self._filter_clonotypes(clonoset, list_type="white")
        if len(self.black_list) > 0:
            clonoset = self._filter_clonotypes(clonoset, list_type="black")
        return clonoset
    
    def _convert_clonoset(self, clonoset, colnames):
        # copy clonoset for not changing the original one
        c_clonoset = clonoset.copy()

        # basic column name in clonoset DF
        
        count_column, fraction_column = decide_count_and_frac_columns(colnames, self.by_umi, suppress_warnings=True)
        
        rename_dict = {count_column: "count",
                       fraction_column: "freq",
                       colnames["cdr3aa_column"]: "cdr3aa",
                       colnames["cdr3nt_column"]: "cdr3nt",
                       colnames["v_column"]: "v",
                       colnames["d_column"]: "d",
                       colnames["j_column"]: "j"}
        c_clonoset = c_clonoset.rename(columns=rename_dict)
        
        result_columns = ["count", "freq"]
        if "cdr3nt" in c_clonoset.columns:
            result_columns.append("cdr3nt")
        result_columns += ["cdr3aa", "v"]
        segment_borders_columns = ["VEnd", "DStart", "DEnd", "JStart"]


        
        
        # In the case of MiXCR and Bioadaptive format the segment type columns
        # usually show several segment variants with particular allele and score.
        # Here we extract only the name of the best hit without allele ane score
        c_clonoset["v"] = c_clonoset["v"].apply(lambda x: extract_segment(x))
        if "d" in c_clonoset.columns:
            c_clonoset["d"] = c_clonoset["d"].apply(lambda x: extract_segment(x))
            result_columns.append("d")
        if "j" in c_clonoset.columns:
            c_clonoset["j"] = c_clonoset["j"].apply(lambda x: extract_segment(x))
            result_columns.append("j")
        
        # add the column for Constant segment if it exists in the original clonoset
        if colnames["c_column"] is not None:
            c_clonoset = c_clonoset.rename(columns={colnames["c_column"]: "c"})
            c_clonoset["c"] = c_clonoset["c"].apply(lambda x: extract_segment(x))
            result_columns += ["c"]
        
        # obtain the borders of the segments within CDR3 region, if possible and add them to
        # resulting clonoset
        if "refPoints" in c_clonoset.columns:
            c_clonoset["VEnd"] = c_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 11, minus=True))
            c_clonoset["DStart"] = c_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 12, minus=False))
            c_clonoset["DEnd"] = c_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 15, minus=True))
            c_clonoset["JStart"] = c_clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 16, minus=False))
        
        result_columns += [col for col in segment_borders_columns if col in c_clonoset.columns]    
        
        # save "sample_id" column if it is present in clonoset
        if "sample_id" in c_clonoset.columns:
            result_columns.append("sample_id")
        c_clonoset = c_clonoset.sort_values(by="count", ascending=False).reset_index(drop=True)
        return c_clonoset[result_columns]
    
    def _make_common_columns(self, clonoset_in, colnames):
        clonoset = clonoset_in.copy()

        # treat refPoints
        refpoints_columns = ["VEnd", "DStart", "DEnd", "JStart"]
        if len(set(clonoset.columns).intersection(set(refpoints_columns))) < 4: # check if not all the columns present
            if "refPoints" in clonoset.columns: # if refPoints is present, create new columns
                clonoset["VEnd"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 11, minus=True))
                clonoset["DStart"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 12, minus=False))
                clonoset["DEnd"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 15, minus=True))
                clonoset["JStart"] = clonoset["refPoints"].apply(lambda x: extract_refpoint_position(x, 16, minus=False))
        
        # treat v,d,j segments

        # V
        if "v" not in clonoset.columns:
            if colnames["v_column"] is not None:
                clonoset["v"] = clonoset[colnames["v_column"]]
        if "v" in clonoset.columns:
            clonoset["v"] = clonoset["v"].apply(lambda x: extract_segment(x))
        
        # D
        if "d" not in clonoset.columns:
            if colnames["d_column"] is not None:
                clonoset["d"] = clonoset[colnames["d_column"]]
        if "d" in clonoset.columns:
            clonoset["d"] = clonoset["d"].apply(lambda x: extract_segment(x))

        # J
        if "j" not in clonoset.columns:
            if colnames["j_column"] is not None:
                clonoset["j"] = clonoset[colnames["j_column"]]
        if "j" in clonoset.columns:
            clonoset["j"] = clonoset["j"].apply(lambda x: extract_segment(x))
        
        # C
        if "c" not in clonoset.columns:
            if colnames["c_column"] is not None:
                clonoset["c"] = clonoset[colnames["c_column"]]
        if "c" in clonoset.columns:
            clonoset["c"] = clonoset["c"].apply(lambda x: extract_segment(x))

        if "cdr3aa" not in clonoset.columns:
            if colnames["cdr3aa_column"] is not None:
                clonoset["cdr3aa"] = clonoset[colnames["cdr3aa_column"]].astype(str)
        if "cdr3nt" not in clonoset.columns:
            if colnames["cdr3nt_column"] is not None:
                clonoset["cdr3nt"] = clonoset[colnames["cdr3nt_column"]].astype(str)

        return clonoset



    def _unweight(self, clonoset_in, colnames):
        clonoset = clonoset_in.copy()
        clonoset[colnames["count_column"]] = 1
        return clonoset
    
    def _recount_fractions_for_clonoset(self, clonoset_in, colnames):
        if self.is_empty():
            return clonoset_in
        clonoset = clonoset_in.copy()

        count_column = colnames["count_column"]
        fraction_column = colnames["fraction_column"]
        umi_column = colnames["umi_column"]
        umi_fraction_column = colnames["umi_fraction_column"]

        clonoset[fraction_column] = clonoset[count_column]/clonoset[count_column].sum()
        if colnames["umi"]:
            clonoset[umi_fraction_column] = clonoset[umi_column]/clonoset[umi_column].sum()
        return clonoset
    
    def _filter_by_count(self, clonoset_in, colnames):
        if self.count_threshold is None:
            return clonoset_in
        clonoset = clonoset_in.copy()
        count_column = colnames["count_column"]
        clonoset = clonoset.loc[clonoset[count_column] >= self.count_threshold]

        return clonoset

    def _check_input(self):

        """
        Check if the object was created properly

        Raises:
            ValueError: in case of incorrect parameter values
        """
        functionality_options = ["a", "f", "n"]
        count_threshold_limits = [0, 100000]
        if self.functionality not in functionality_options:
            raise ValueError(f"Incorrect value '{self.functionality}' for functionality. Possible values: {', '.join(functionality_options)}")
        if self.count_threshold is not None:
            if not isinstance(self.count_threshold, int):
                raise TypeError("Count threshold must be an 'int' or 'None'")
            if (self.count_threshold < count_threshold_limits[0]
                  or self.count_threshold > count_threshold_limits[1]):
                raise ValueError(f"Incorrect value '{self.functionality}' for count_threshold. Possible values: {count_threshold_limits}")                
        if self.downsample_size is not None:
            if not isinstance(self.downsample_size, int):
                raise ValueError(f"Incorrect value '{self.downsample_size}' for downsample_size. Only int or None possible")
            elif self.downsample_size < 1:
                raise ValueError(f"Incorrect value '{self.downsample_size}' for downsample_size. Value too low")
        if self.top is not None:
            if not isinstance(self.top, int):
                raise ValueError(f"Incorrect value '{self.top}' for top. Only int or None possible")
            elif self.top < 1:
                raise ValueError(f"Incorrect value '{self.top}' for top. Value too low")
        if not isinstance(self.seed, Hashable):
            raise ValueError(f"Incorrect value '{self.seed}' for seed. Must be hashable")
        pool_by_options = ["", "aa", "aaV", "aaVj", "nt", "ntV", "ntVJ"]
        if self.pool_by not in pool_by_options:
            raise ValueError(f"Incorrect value '{self.pool_by}' for clonoset pool. Possible values: {', '.join(pool_by_options)}")
  
    def _downsample(self, clonoset_in, colnames):
        """
        Downsample clonoset.
        
        This function takes the total number of reads or UMIs of the clonoset.
        Then randomly samples the downsample_size from 0 to this total number of reads/UMIs. 
        This random sample is mapped to the clonotype sizes and 
        the new downsampled clonoset is created
        """
        
        if self.downsample_size is None:
            return clonoset_in
        
        clonoset = clonoset_in.copy()
        count_column = colnames["count_column"]
        total_count = int(clonoset[count_column].sum())

        # raise ValueError if UMI/read count is less then downsample_size
        
        if total_count < self.downsample_size:
            if self.ignore_small_clonosets:
                return clonoset
            else:
                raise ValueError(f"total count {total_count} is less than downsample size {self.downsample_size}")
        elif total_count == self.downsample_size:
            return clonoset
        
        # set seed if given and take the sample of total_count
        if self.seed is not None:
            random.seed(self.seed)
        sample = sorted(random.sample(range(total_count), self.downsample_size))
        
        # map the sample to the clone counts in the clonoset
        curr_sum = 0
        i = 0
        new_counts_dict = {}
        for index,r in clonoset.iterrows():
            curr_sum+=r[count_column]
            new_count = 0
            if i == self.downsample_size:
                break
            while(sample[i]<curr_sum):
                new_count+=1
                i+=1
                if i == self.downsample_size:
                    break
            if new_count > 0:
                new_counts_dict[index]=new_count
        
        # filter clonoset for missed clones and set new clone counts
        (indices,counts) = zip(*new_counts_dict.items())
        clonoset = clonoset.loc[clonoset.index.isin(indices)]
        clonoset[count_column] = counts    
        return clonoset.reset_index(drop=True)
    
    def _get_top(self, clonoset_in, colnames):
        """
        Takes top N biggest clones from the clonoset.
        
        Mix-tails is recommended for use, because the order of the clonotypes
        with equal count may not be independent from their other properties.
        This option mixes up the order of all clonotypes in clonoset and then
        sorts them by count in decreasing order, so that clonotypes with the same
        count not have completely random order. Also use seed option for reproducibility
        of the results.

        Raises:
            ValueError: if clone count is less then required top
        """

        if self.top is None:
            return clonoset_in
        
        clonoset = clonoset_in.copy()
        count_column = colnames["count_column"]
        
        #shuffle the order of clonotypes if required
        if self.seed is not None:
            random.seed(self.seed)
        if self.mix_tails:
            index_order = random.sample(list(clonoset.index), len(clonoset))
            clonoset = clonoset.iloc[index_order] 
            clonoset = clonoset.sort_values(by=count_column, ascending=False)

        if self.top > len(clonoset):
            if self.ignore_small_clonosets:
                return clonoset
            else:
                raise ValueError(f"Warning! Clonoset size - {len(clonoset)} - is less than required top - {self.top}")
        
        # take top
        if self.top > 0:
            clonoset=clonoset.iloc[:self.top]

        return clonoset.reset_index(drop=True)
    
    def _filter_by_functionality(self, clonoset_in, colnames):
        clonoset = clonoset_in.copy()
        
        if self.functionality == "f":
            clonoset = filter_by_functionality(clonoset, colnames=colnames, functional=True)
        if self.functionality == "n":
            clonoset = filter_by_functionality(clonoset, colnames=colnames, functional=False)
            
        return clonoset.reset_index(drop=True)
            
    def __str__(self):
        
        functionality = {"a": "any",
                         "f": "only functional clones (no frameshifts and Stops)",
                         "n": "only non-functional"}
        output  = f"Filter name:\t{self.name}\n"
        output += f"Functionality:\t{functionality[self.functionality]}\n"
        random = False
        change_size = False
        if self.functionality != "a":
            change_size = True
        if isinstance(self.count_threshold, int):
            output += f"Count threshold:\t{self.count_threshold}\n"
            change_size = True
        if isinstance(self.downsample_size, int):
            output += f"Downsample size:\t{self.downsample_size}\n"
            random = True
            change_size = True
        if isinstance(self.top, int):
            output += f"Take top:\t{self.top}\n"
            change_size = True
            if self.mix_tails:
                random = True
        if self.unweight:
            output += f"Clone size unweighted (all clone counts = 1)"
                
        if change_size:
            if self.by_umi:
                output += f"Count by:\t UMI (if exist)\n"
            else:
                output += f"Count by:\t reads/counts\n"
        if random:
            if isinstance(self.seed, int):
                output += f"Seed for random:\t{self.seed}\n"
            else:
                output += f"Seed for random:\tunset\n"
                output += f"Warning: filter contains random events.\nTo obtain reproducible results, set seed in calcutations or manually (random.seed(some_int))\nprior to applying the filter.\nNote only seed that is set as this object parameter will work for mix_tails\n"
                 
        return output
    
    def _pool_clonoset(self, clonoset_in, colnames):
        # copy clonoset and sort by clone counts and reset index for order
        clonoset = clonoset_in.copy().sort_values(by=colnames["count_column"], ascending=False).reset_index(drop=True)
        
        # create list of pool columns
        aa, check_v, check_j = overlap_type_to_flags(self.pool_by)
        columns_for_pool = []
        if aa:
            columns_for_pool.append(colnames["cdr3aa_column"])
        else:
            columns_for_pool.append(colnames["cdr3nt_column"])
        if check_v:
            columns_for_pool.append(colnames["v_column"])
        if check_j:
            columns_for_pool.append(colnames["j_column"])

        # create column combining all pool columns
        clonoset["pool_id"] = clonoset.apply(lambda x: "|".join([x[colname] for colname in columns_for_pool]), axis=1)
        
        indices_to_retain = []

        for pool_id in clonoset["pool_id"].unique():
            pool_clonoset = clonoset.loc[clonoset["pool_id"] == pool_id]

            # select the clone with biggest count - it will represent pooled clonotypes by
            # columns other that count and freq
            top_index = pool_clonoset.index[0]
            indices_to_retain.append(top_index)
            
            # sum counts and fractions for pooled clonotypes
            clonoset.loc[top_index,colnames["count_column"]] = pool_clonoset[colnames["count_column"]].sum()
            clonoset.loc[top_index,colnames["fraction_column"]] = pool_clonoset[colnames["fraction_column"]].sum()
        
        # retain only rows with representative clonotypes and remove technical column
        clonoset = clonoset.loc[indices_to_retain].drop(columns=["pool_id"])

        return clonoset

    def _filter_clonotypes(self, clonoset_in, list_type):
        if list_type == "white":
            clonotypes_list = self.white_list
        elif list_type == "black":
            clonotypes_list = self.black_list
        else:
            raise ValueError("list_type must be 'white' or 'black'")
        
        clonoset = clonoset_in.copy()
        
        clonoset["filter_pass"] = clonoset.apply(lambda x: self._compare_clonoset_list_row_with_clonotype(x, clonotypes_list), axis=1)
        if list_type == "white":
            clonoset = clonoset[clonoset["filter_pass"]].drop(columns=["filter_pass"]).reset_index(drop=True)
        else:
            clonoset = clonoset[~clonoset["filter_pass"]].drop(columns=["filter_pass"]).reset_index(drop=True)
        return clonoset


# def convert_clonoset_to_clonotype_filter_list(clonoset_df, overlap_type="aaVJ"):
#     clonotypes_list = []
#     aa, include_v, include_j = intersections.overlap_type_to_flags(overlap_type)
#     for i,r in clonoset_df.iterrows():
#         clonotype = []
#         if aa:
#             clonotype.append(row["cdr3aa"])
#         else:
#             clonotype.append(row["cdr3nt"])
#         if include_v:
#             clonotype.append(row["v"])
#         if include_j:
#             clonotype.append(row["j"])
    
#         clonotype = tuple(clonotype)
#         clonotypes_list.append(clonotype)
#     return clonotypes_list
     
    def _compare_clonoset_row_with_clonotype(self, row, clonotype):
        c_len = len(clonotype)
        if c_len == 1:
            if row["cdr3aa"] == clonotype[0]:
                return True
        elif c_len == 2:
            if row["cdr3aa"] == clonotype[0] and row["v"] == clonotype[1]:
                return True
        elif c_len == 3:
            if row["cdr3aa"] == clonotype[0] and row["v"] == clonotype[1] and row["j"] == clonotype[2]:
                return True
        else:
            # need to write better explanation for error
            raise ValueError("clonotypes must contain from 1 to 3 values")
            
        return False

    def _compare_clonoset_list_row_with_clonotype(self, row, clonotypes_list):
        for clonotype in clonotypes_list:
            if self._compare_clonoset_row_with_clonotype(row, clonotype):
                return True
        return False

    def is_empty(self):
        return self.functionality == "a" and self.downsample_size is None and self.top is None and self.count_threshold is None and not self.unweight

    def _repr_html_(self):
        """
        function for printing the Filter properties to Jupyter output
        """

        functionality = {"a": "any",
                         "f": "only functional clones (no frameshifts and Stops)",
                         "n": "only non-functional"}
        output  = f"<p>Filter name: {self.name}</p>"
        output += f"<p>Functionality:\t{functionality[self.functionality]}</p>"
        random = False
        change_size = False
        if self.functionality != "a":
            change_size = True
        if isinstance(self.count_threshold, int):
            output += f"<p>Count threshold: {self.count_threshold}</p>"
            change_size = True
        if isinstance(self.downsample_size, int):
            output += f"<p>Downsample size: {self.downsample_size}</p>"
            random = True
            change_size = True
        if isinstance(self.top, int):
            output += f"<p>Take top: {self.top}</p>"
            change_size = True
            if self.mix_tails:
                random = True
        if self.unweight:
            output += f"<p>Clone size unweighted: (all clone counts = 1)</p>"
                
        if change_size:
            if self.by_umi:
                output += f"<p>Count by:  UMI (if exist)</p>"
            else:
                output += f"<p>Count by: reads/counts</p>"
        if random:
            if isinstance(self.seed, int):
                output += f"<p>Seed for random: {self.seed}</p>"
            else:
                output += f"<p>Seed for random:\tunset</p>"
                output += f"<p>Warning: filter contains random events. To obtain reproducible results, set seed in calcutations or manually (random.seed(some_int)) prior to applying the filter. Note only seed that is set as this object parameter will work for mix_tails</p>"
                 
        return output
