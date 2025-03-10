---
layout: post
title: "Report from Geomorphometry 2009"
date: "2009-11-01"
tags: [story,event,Zurich,Switzerland]
hide_hero: true
published: true
image: false
author: Bob MacMillan
---

FYI: a short report from Geomorphometry 2009 by Bob MacMillan (ISRIC) published in the Pedometron newsletter #28:

"The main purpose of my participation was to keep informed about other efforts, similar to GlobalSoilMap.net, that have a interest in processing digital elevation data and other digital data sets globally or at least for extremely large areas. This conference actually contained a large number of presentations of direct relevance for the GlobalSoilMap.net project. Perhaps first and foremost were the descriptions of efforts being undertaken in Australia (Gallant and Read) and Europe (Köthe and Bock) to process SRTM DEM data at 30 m (Australia) and 90 m (Europe) grid resolution to reduce artefacts and produce a filtered and cleaned DEM that is more suitable for use to produce inputs for the GlobalSoilMap.net project. Both of these presentations highlighted the significant advantages that can be realised by applying a series of filtering and conditioning routines to the original raw SRTM DEM data. It is obvious that similar procedures would prove equally useful if applied to SRTM DEM data sets for other parts of the world under the jurisdiction of other GlobalSoilMap.net nodes. Gallant has offered to help with efforts in other Nodes if asked.

Also of great interest were several projects that demonstrated that it is indeed possible to process and produce digital output for global scale digital data sets, including global scale SRTM DEM data sets. Reuter and Nelson presented a description of WorldTerrain, a contribution of the Global Geomorphometric Atlas. Peter Guth described processing of global scale SRTM data to identify and classify organized linear landforms (dunes). Peter also provided examples of multiple scale analysis and illustrated what you get to “see” from DEMs of 1 m, 100 m and 2 km grid resolution. Guth intends to publish the many different grids of DEM derivatives he produced for his project and make these processed data available for free and widespread use by others. Marcello Gorini described a physiographic classification of the ocean flood using a multi-resolution geomorphometric approach.

Several authors presented methods that may prove of interest to the GlobalSoilMap.net project. Gallant and Hutchinson described a differential equation for computing specific catchment area that could be applied to produce an improved terrain covariate for use in the GlobalSoilMap.net project. Similarly, Peckham, gave a new algorithm for creating DEMs with smooth elevation profiles that could be used to condition rough SRTM or GDEM data sets to smooth out noise and produce more hydrologically plausible surfaces. This algorithm was of particular interest to the GlobalSoilMap.net project because it appeared to be able to introduce hydrologically and geomorphologically relevant detail into 90 m SRTM DEMs of relatively low spatial detail.

Romstad and Etzelmuller described a new approach for segmenting hillslopes into landform elements by applying a watershed algorithm to a surface defined by the total curvature at a point instead of the raw elevation value. The resulting watersheds were bounded by lines of maximum curvature, effectively structuring each hillslope into components partitioned by lines of maximum local curvature. This is harder to explain than to understand when illustrated but it is remarkably simple to implement and may provide a new way of automatically segmenting hillslopes in a simple and efficient fashion.

Metz and others presented an algorithm for fast and efficient processing of massive DEMs to extract drainage networks and flow paths. This is of considerable interest and relevance to the GlobalSoilMap.net project because of the project’s need to process SRTM data globally to compute hydrological flow networks and various indices that are computed based on flow networks (e.g. elevation above channel, distance from divide). This algorithm can process data sets of hundreds of millions of cells (11,424 rows by 13,691 cols) in a few minutes instead of a few days (or not at all for some algorithms that fail on data sets this large).

Overall, this was an excellent conference, dominated by leading edge research in the area of geomorphic processing of digital elevation data that is of direct relevance and interest to the GlobalSoilMap.net project. We have much to learn from these researchers and much to benefit from maintaining contacts and working relationships with them."
