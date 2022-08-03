This folder contains the parameter files used for benchmarking Comet and MSFragger against Carina

I tried to make the comparison as close as possible (generally using default parameters where applicable), but there are a couple caveats:

- Comet doesn't have an easily configuration fragment PPM tolerance, instead it uses a bin width (this is used for integerizing all m/z values)
- I turned off mass calibration and parameter optimization in MSFragger, since it can take quite a long time and had minimal benefits on the high-quality datasets I used