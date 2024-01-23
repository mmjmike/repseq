from .clonosets import (decide_count_and_frac_columns,
                        get_column_names_from_clonoset,
                        filter_by_functionality)

from .common_functions import (extract_segment,
                               extract_refpoint_position)

from collections import Hashable

import random

class Filter:
    
    """
    Clonoset filter.
    May be used to filter clonosets:
        - by clone functionality
        - randomly downsample them to particular number of reads of UMIs
        - take top clonotypes by size (number of reads of UMIs) with or without random mixing

    Args:
        - name (str): the name of the filter. Will be displayed in print
        - functionality (str): Possible values:
            - "a" - any (default). No clones are filtered out
            - "f" - only functional. Those, not having stop codons and 
                frame shifts in CDR3 regions, or having non-empty values in CDR3 amino-acid
                sequence
            - "n" - only-nonfunctional - opposite to "f" - functional
        - downsample (int): the number of reads/UMIs to randomly downsample the clonoset to.
            default value 'None' - means not to apply downsampling
        - top (int): the number of top biggest by reads/UMIs clonotypes to take from the clonoset.
            default value 'None' - means not to apply top
        - by_umi (bool): default=False. Which column to take for clonotype count - reads or UMIs 
            (if UMI count column exists).
        - mix_tails (bool): default=False. Defines whether to randomly mix-up the order of clonotypes
            before sorting by size and taking the top clonotypes. Basically mix_tails=True mixes up 
            clonotypes with the same size in read or UMIs.
        - count_threshold (int): limits [0:100000], all clonotypes with count less than this value will
            be filtered out
        - seed (any hashable type): better to use int - seed for reproducibility of random events 
            (downsampling or top with mix-tails). Default=None.
    """

    def __init__(self, name="default_filter", functionality="a", downsample=None,
                 top=None, by_umi=False, mix_tails=False, count_threshold=None, seed=None):
        self.name = name
        self.functionality = functionality
        self.downsample_size = downsample
        self.top = top
        self.by_umi = by_umi
        self.mix_tails = mix_tails
        self.seed = seed
        self.count_threshold = count_threshold
        self._check_input()
        
    def spawn(self):
        """
        
        Returns:
            the copy of the filter. Necessary for parallel computing

        """
        return Filter(name=self.name, functionality=self.functionality,
                      downsample=self.downsample_size, top=self.top,
                      by_umi=self.by_umi, mix_tails=self.mix_tails,
                      count_threshold=self.count_threshold, seed=self.seed)
        
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

        # converting to common VDJtools-like format and obtaining new colnames
        clonoset = self._convert_clonoset(clonoset, colnames)
        colnames = get_column_names_from_clonoset(clonoset)

        # application of main filters
        if self.functionality != "a":
            clonoset = self._filter_by_functionality(clonoset, colnames)
        
        clonoset = self._filter_by_count(clonoset, colnames)
        
        clonoset = self._downsample(clonoset, colnames)
        clonoset = self._get_top(clonoset, colnames)

        # the fraction columns need to be recounted after filtering, as they
        # remain the same as in the original clonoset before filtration
        clonoset = self._recount_fractions_for_clonoset(clonoset, colnames)
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
        
        return c_clonoset[result_columns]
        
    
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
        if functionality != "a":
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
    
    def is_empty(self):
        return self.functionality == "a" and self.downsample_size is None and self.top is None and self.count_threshold is None

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
        if functionality != "a":
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
