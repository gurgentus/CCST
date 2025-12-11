# FastAPI router for data operations
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
import data  # Import the pure Python data module

router = APIRouter()

# Request models
class DataRequest(BaseModel):
    file: str
    name: str

class FitRequest(BaseModel):
    name: str
    predictor: str

class GaussianProcessRequest(BaseModel):
    data_name: str
    plot: str

class UniformGPRequest(BaseModel):
    plot: str

class PredictRequest(BaseModel):
    data: float
    predictor: str

class PlotDataRequest(BaseModel):
    name: str
    plot: str


# Endpoints - thin layer that calls data module functions
@router.post("/data")
async def data_endpoint(req: DataRequest):
    """Associate a file with a variable name"""
    try:
        result = data.data(req.file, req.name)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/fit")
async def fit_endpoint(req: FitRequest):
    """Fit a linear regression model to data"""
    try:
        data.fit(req.name, req.predictor)
        return {"result": "Linear regression model fitted"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/gaussian_process")
async def gaussian_process_endpoint(req: GaussianProcessRequest):
    """Generate gaussian process plot from data"""
    try:
        result = data.gaussian_process(req.data_name, req.plot)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/uniform_gaussian_process")
async def uniform_gp_endpoint(req: UniformGPRequest):
    """Generate uniform gaussian process plot"""
    try:
        result = data.uniform_gaussian_process(req.plot)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/predict")
async def predict_endpoint(req: PredictRequest):
    """Make a prediction using a trained predictor"""
    try:
        result = data.predict(req.data, req.predictor)
        return {"result": float(result)}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/plotdata")
async def plotdata_endpoint(req: PlotDataRequest):
    """Plot data from a CSV file"""
    try:
        result = data.plotdata(req.name, req.plot)
        return {"result": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
