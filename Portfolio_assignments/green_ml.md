# Assignment: Assessing the Environmental Impact of Your ML Application

The rapid growth and adoption of machine learning (ML), particularly deep neural networks (DNNs), have significantly influenced the evolution of computer architectures. The demand for high-performance computations has driven innovations in hardware, including the development of specialized processors like Tensor Processing Units (TPUs), Graphics Processing Units (GPUs), and advancements in multi-core Central Processing Units (CPUs). These technologies aim to meet the computational needs of increasingly complex ML models.

To optimize performance, modern solutions leverage efficient linear algebra kernels provided by frameworks like TensorFlow, NVIDIA's CUDA-X, and other specialized libraries. These tools enable rapid matrix computations, the backbone of DNN training and inference. Additionally, trade-offs in numerical precision have become a key strategy for balancing accuracy and efficiency. While higher precision using more bits can improve model fidelity, using fewer bits, such as integers, can reduce memory consumption and power usage without significant loss in performance.

As ML applications continue to scale, understanding and mitigating their environmental impact becomes increasingly important. This assignment explores how these architectural innovations, while boosting computational efficiency, influence the energy consumption and carbon footprint of ML applications. Students will evaluate their own ML projects to assess sustainability and propose strategies to optimize both performance and environmental impact.

## Objective
Evaluate the environmental impact of your ML or AI application and propose strategies to minimize its carbon footprint. Use either one of your existing ML applications or your current Omics project application. 

## Preparatory Reading:
Read (as a minimum) the following two articles (more articles can be found at the bottom of this page): 
L. Lannelongue, J. Grealey, M. Inouye, Green Algorithms: Quantifying the Carbon Footprint of Computation. Adv. Sci. 2021, 8, 2100707. https://doi.org/10.1002/advs.202100707

Lannelongue L, Grealey J, Bateman A, Inouye M (2021) Ten simple rules to make your computing more environmentally sustainable. PLoS Comput Biol 17(9): e1009324. https://doi. org/10.1371/journal.pcbi.1009324

## Assignment
Give a brief description of your ML application, including: 
-	Model type
-	Dataset size
-	Training setup (compute type, hardware used)
-	Deployment setup (if applicable)
Calculate the approximate energy consumption of your ML application using tools like: 
-	https://calculator.green-algorithms.org/
-	CodeCarbon
-	MLCO2 Impact Calculator

Suggest optimizations to make your ML application more sustainable

## Deliverables: include in your portfolio A 2-3 page report with your steps findings: 
-	Model description
-	Results of your environmental impact assessment.
-	Strategies for optimizing the application’s sustainability.


## references
Han, Song, Xingyu Liu, Huizi Mao, Jing Pu, Ardavan Pedram, Mark A. Horowitz, and William J. Dally. “EIE: Efficient Inference Engine on Compressed Deep Neural Network.” SIGARCH Comput. Archit. News 44, no. 3 (June 2016): 243–54. https://doi.org/10.1145/3007787.3001163.

Horowitz, Mark. “1.1 Computing’s Energy Problem (and What We Can Do about It).” In 2014 IEEE International Solid-State Circuits Conference Digest of Technical Papers (ISSCC), 10–14. IEEE, 2014.

Lannelongue, Loïc, Jason Grealey, and Michael Inouye. “Green Algorithms: Quantifying the Carbon Footprint of Computation.” Advanced Science 8, no. 12 (2021): 2100707.

Lannelongue L, Grealey J, Bateman A, Inouye M (2021) Ten simple rules to make your computing more environmentally sustainable. PLoS Comput Biol 17(9): e1009324. https://doi. org/10.1371/journal.pcbi.1009324 

Martinez, D.R. and Kifle, B.M., 2024. Artificial Intelligence: A Systems Approach from Architecture Principles to Deployment. MIT Press.

Patterson, David A., Joseph Gonzalez, Quoc V. Le, Chen Liang, Lluis-Miquel Munguia, Daniel Rothchild, David R. So, Maud Texier, and Jeff Dean. “Carbon Emissions and Large Neural Network Training.” CoRR abs/2104.10350 (2021). https://arxiv.org/abs/2104.10350.

Sze, Vivienne, Yu-Hsin Chen, Tien-Ju Yang, and Joel S. Emer. “Efficient Processing of Deep Neural Networks: A Tutorial and Survey.” Proceedings of the IEEE 105, no. 12 (2017): 2295–2329. https://doi.org/10.1109/JPROC.2017.2761740.
