### Fine-grained parallelization approach

**TL;DR: The first parallelization approach (upfront work splitting) was very coarse-grained and could cause some nodes to take much longer than others, resulting in timeouts on synchronization. Therefore it was switched to a as fine-grained as possible approach.**

Assume **several workers** and one master.

We can run a batch of **N photons in parallel**, where N equals thread count. In order to limit time spent on GPU, we always run for a fixed **number of bounces**. We cannot know how many runs it will take until the batch is finished (all photons terminated). Instead we count finished photons after each run. If after run **i**, a number **K > Z photons have finished**, the worker sends K to master. One thread of the master must be in receiving state. If the master finds that the **number of photons that is not being worked on is L >= K**, the response is **CONTINUE**, otherwise **FINISH**. When receiving CONTINUE, the worker spawns K new photons, otherwise it finishes the remaining photons (with thread utilitization < N) and calls **Reduce** on its output arrays. When L drops to zero, the master also calls **Reduce**, collects the total result and writes it to file.

Important: After the master responds with FINISH, he should finish the remaining L photons on his own, because there is no guarantee that a worker has space for exactly L photons at some point, meaning (in the negative case) they could either run too many or L never drops to zero.

Note: It is expected that Reduce doesn't cause timeouts, since all workers receive a FINISH response within a short amount of time, because of the high communication frequency.

The high communication frequency leads to a tradeoff between thread utilitization and communication efficiency, both impacting performance to unknown extent. If we **set Z to zero**, we get maximum thread utilitization, but high communication overhead. The **higher Z** becomes, the worse thread utilitization and the better communication efficiency becomes.