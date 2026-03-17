
import narwhals
import scipy.stats
from attrs import frozen

from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import PhenotypeInfo, QuantPhenotype
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import StandardMagmaTaskGenerator


def _get_liability_scale_heritability(
    observed_heritability: float,
    sample_prev: float,
    population_prev: float,
):
    T= scipy.stats.norm.ppf(1-population_prev)
    a= scipy.stats.norm.pdf(T)
    return observed_heritability*(population_prev*(1-population_prev))**2/a**2/(sample_prev*(1-sample_prev))


@frozen
class HeritabilityConversionPipe(DataProcessingPipe):
    observed_heritability_column:str
    liability_heritability_column:str
    sample_info : PhenotypeInfo
    assume_even_sample: bool # useful if we are working with an effective sample size, and so have implicitly converted to sample-prev=0.5

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        if isinstance(self.sample_info, QuantPhenotype ):
            return x
        df = x.collect().to_pandas()
        if self.assume_even_sample:
            sample_prev=0.5
        else:
            sample_prev =self.sample_info.sample_prevalence
        df[self.liability_heritability_column] = df[self.observed_heritability_column].apply(
        lambda h:
            _get_liability_scale_heritability(h, sample_prev=sample_prev, population_prev=self.sample_info.estimated_population_prevalence)
        )
        return narwhals.from_native(df).lazy()