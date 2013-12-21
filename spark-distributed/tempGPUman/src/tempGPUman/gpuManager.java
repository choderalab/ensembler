package tempGPUman;
import java.util.concurrent.locks.ReentrantLock;
import py4j.GatewayServer;



public class gpuManager {

	private int numGPUs, workerPerGPU;
	public int[] GPUlist;
	//we're doing to use this synchronization primitive, since it guarantees "fairness"
	//otherwise, I'm afraid that the relinquishers may never get a chance to do so
	private final ReentrantLock lock = new ReentrantLock(true);

	public gpuManager(int numGPU, int workersPerGPU)
	{
		numGPUs = numGPU;
		workerPerGPU = workersPerGPU;
		GPUlist = new int[numGPUs];
		
	}
	
	//this is called in order to get an assigned GPU number
	public int getGPUAssignment()

	
	{
		
		boolean assigned = false;
		//make sure that two people are not looking/writing at the same time;
		//hopefully this works
		while (!assigned){
			
		lock.lock();	
		try {
			for (int i=0; i<numGPUs;i++)
			{
				if (GPUlist[i]<workerPerGPU)
				{
					GPUlist[i]++;
					return i;
				}
				
			}
		}
		finally{
			lock.unlock();
		}
		
			
		}
		//should never reach here; returning -1 means that a GPU was not assigned
		return -1;
		
	}
	
	
	//this is called when the caller is finished with the GPU
	//there should probably be a better way of doing this than relying on
	//the caller to be scrupulous enough to 
	public void relinquishGPUAssignment(int GPUID)
	{
		lock.lock();
		try
		{
			GPUlist[GPUID]--;
		}
		finally{
			lock.unlock();
		}
		
	}
	

	
	
	public static void main(String[] args) {
		
		//start the gateway server, use port 20000 because the spark gateway will also be running
		//this is not a permanent solution to GPU management, but will serve as a stopgap measure
		
		
		//here we take as default 4 GPUs, 1 worker/GPU (common configuration on the cluster)
		//I will add arguments to this
	
		int numGPUs =4;
		int workersPerGPU =1;
	
		System.out.println("Starting gateway server with "+ numGPUs + " gpus and " + workersPerGPU + " workers/gpu");
	    GatewayServer gatewayServer = new GatewayServer(new GPUManagerEntry(numGPUs,workersPerGPU),20000);
	    gatewayServer.start();
	    System.out.println("Started gateway server");

		
		

	}

}
